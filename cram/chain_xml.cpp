#include "cram/chain_xml.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#ifdef CRAM_WITH_CHAIN_XML
#include <pugixml.hpp>
#endif

namespace cram {
namespace {

// Element symbols indexed by proton number; index 0 is a placeholder.
constexpr std::array<std::string_view, 119> kSymbols = {
    "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf",
    "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

int symbolToZ(std::string_view sym) {
  for (int z = 1; z < static_cast<int>(kSymbols.size()); ++z)
    if (kSymbols[static_cast<std::size_t>(z)] == sym)
      return z;
  return -1;
}

}  // namespace

std::string elementSymbol(int z) {
  if (z < 1 || z >= static_cast<int>(kSymbols.size()))
    return {};
  return std::string(kSymbols[static_cast<std::size_t>(z)]);
}

std::optional<Zai> parseNuclideName(std::string_view name) {
  std::size_t i = 0;
  while (i < name.size() && (std::isalpha(static_cast<unsigned char>(name[i])) != 0))
    ++i;
  if (i == 0 || i >= name.size())
    return std::nullopt;
  const int z = symbolToZ(name.substr(0, i));
  if (z < 0)
    return std::nullopt;

  // Mass number digits follow the symbol.
  std::size_t j = i;
  while (j < name.size() && (std::isdigit(static_cast<unsigned char>(name[j])) != 0))
    ++j;
  if (j == i)
    return std::nullopt;
  int a = 0;
  std::from_chars(name.data() + i, name.data() + j, a);

  // Optional metastable suffix "_m<level>".
  int iso = 0;
  if (j < name.size() && name[j] == '_') {
    std::size_t k = j + 1;
    if (k < name.size() && (name[k] == 'm' || name[k] == 'M'))
      ++k;
    if (k < name.size())
      std::from_chars(name.data() + k, name.data() + name.size(), iso);
  }

  return Zai{.z = z, .a = a, .i = iso};
}

#ifdef CRAM_WITH_CHAIN_XML

namespace {

// Whitespace-tokenize element text into a vector of strings.
std::vector<std::string> tokenize(const char* text) {
  std::vector<std::string> out;
  if (text == nullptr)
    return out;
  std::istringstream is(text);
  std::string tok;
  while (is >> tok)
    out.push_back(tok);
  return out;
}

void parseDecay(const pugi::xml_node& nuc, const Zai& parent, DepletionChain& chain,
                ChainXmlDiagnostics& diag) {
  const double halfLife = nuc.attribute("half_life").as_double(0.0);
  if (halfLife <= 0.0)
    return;  // stable (or no half-life given)

  DecayData d{.halfLife = halfLife};
  for (pugi::xml_node dn : nuc.children("decay")) {
    const double branching = dn.attribute("branching_ratio").as_double(0.0);
    const std::string_view type = dn.attribute("type").as_string();
    const std::string_view target = dn.attribute("target").as_string();

    if (target.empty() || target == "Nothing") {
      // OpenMC omits the target whenever the daughter falls outside the chain.
      // The mode is not built and its branching is lost, but the matrix is
      // still right: the parent is removed at the full decay constant either
      // way, and nothing is produced. Not a parse failure -- parseReactions()
      // reads the same two spellings the same way -- so it is not counted.
      //
      // A targetless "sf" is different: it is this library's own spelling for
      // spontaneous fission with products from the SFY table (yields at energy
      // 0). rtyp 6.0 is the ENDF code for it; decayDaughter() never consults
      // rtyp for a fission mode, but it keeps the mode self-describing.
      if (type == "sf")
        d.modes.push_back(
            DecayMode{.rtyp = 6.0, .branching = branching, .finalState = 0, .isFission = true});
      continue;
    }

    const std::optional<Zai> daughter = parseNuclideName(target);
    if (!daughter) {
      ++diag.unparsedDecayTargets;  // mode dropped; its branching is lost
      continue;
    }
    // Every named target routes its branch the same way, `sf` included: an
    // OpenMC chain names the parent itself there, since the file does not
    // track SF products as such, and its matrix carries that target like any
    // other daughter, so the branch removes and restores the same atoms.
    // Routing it to the yield tables instead would remove the parent at
    // branching*lambda and substitute the nearest NEUTRON-induced yield set as
    // the products -- and where the parent has no table at all, the atoms
    // would simply vanish.
    //
    // rtyp 0 is deliberate: the explicit daughter is what routes the branch,
    // and matrix assembly never derives one from rtyp when a daughter is set.
    d.modes.push_back(DecayMode{.rtyp = 0.0,
                                .branching = branching,
                                .finalState = daughter->i,
                                .isFission = false,
                                .daughter = daughter});
  }
  chain.setDecay(parent, std::move(d));
}

// Store the tables of one <neutron_fission_yields> block on `parent`. The
// block is not necessarily `parent`'s own: a delegating nuclide is given the
// tables of the nuclide it delegates to (see resolveYieldDelegation()).
void parseFissionYields(const pugi::xml_node& nfy, const Zai& parent, DepletionChain& chain,
                        ChainXmlDiagnostics& diag) {
  for (pugi::xml_node fy : nfy.children("fission_yields")) {
    FissionYields y{.energy = fy.attribute("energy").as_double(0.0)};
    const auto products = tokenize(fy.child("products").text().get());
    const auto data = tokenize(fy.child("data").text().get());
    const std::size_t n = std::min(products.size(), data.size());
    for (std::size_t k = 0; k < n; ++k) {
      const std::optional<Zai> prod = parseNuclideName(products[k]);
      if (!prod) {
        ++diag.unparsedYieldProducts;
        continue;
      }
      y.products.emplace_back(*prod, std::stod(data[k]));
    }
    chain.addFissionYields(parent, std::move(y));
  }
}

// The <neutron_fission_yields> block holding the actual tables for the nuclide
// named `name`, following a chain of parent="..." delegations. Null when the
// name is not in the file, when it has no yield block, or when the delegations
// do not terminate (a cycle, or a longer chain than any real file has).
pugi::xml_node resolveYieldDelegation(const std::unordered_map<std::string, pugi::xml_node>& byName,
                                      std::string name) {
  constexpr int kMaxHops = 8;
  for (int hop = 0; hop < kMaxHops; ++hop) {
    const auto it = byName.find(name);
    if (it == byName.end())
      return {};
    pugi::xml_node nfy = it->second.child("neutron_fission_yields");
    if (!nfy)
      return {};
    const std::string_view next = nfy.attribute("parent").as_string();
    if (next.empty())
      return nfy;  // the tables themselves
    name = next;
  }
  return {};
}

void parseReactions(const pugi::xml_node& nuc, const Zai& parent, DepletionChain& chain,
                    std::vector<ChainReaction>& out, ChainXmlDiagnostics& diag) {
  for (pugi::xml_node rn : nuc.children("reaction")) {
    const std::optional<ReactionType> type = reactionTypeFromName(rn.attribute("type").as_string());
    if (!type) {
      ++diag.unmodeledReactions;
      continue;
    }
    const std::string_view target = rn.attribute("target").as_string();
    std::optional<Zai> product;
    if (!target.empty() && target != "Nothing") {
      product = parseNuclideName(target);
      if (product)
        // Register it: a file may name a product it never declares as its own
        // <nuclide> element, and an unregistered target is one assemble()
        // consumes the parent for and produces nothing from -- the same
        // treatment as a deliberate "Nothing", with nothing to tell them
        // apart. This is what close() does for an unregistered decay daughter,
        // and add() is idempotent, so a declared target costs a hash lookup.
        chain.add(*product);
      else
        ++diag.unparsedReactionTargets;  // product lost; the parent is still consumed
    }
    out.push_back(ChainReaction{.parent = parent,
                                .type = *type,
                                .target = product,
                                .q = rn.attribute("Q").as_double(0.0),
                                .branching = rn.attribute("branching_ratio").as_double(1.0)});
  }
}

// Report on stderr what the reader could not use, one line per condition that
// actually occurred. Only called when the caller passed no diagnostics pointer.
void warnDiagnostics(const ChainXmlDiagnostics& diag, const std::string& what) {
  if (diag.clean())
    return;
  const std::pair<int, const char*> counts[] = {
      {diag.unparsedNuclides, "nuclide name(s) that did not parse (the whole element skipped)"},
      {diag.unparsedDecayTargets, "decay target(s) that did not parse (the mode skipped)"},
      {diag.unparsedYieldProducts, "fission-yield product(s) that did not parse"},
      {diag.unmodeledReactions, "reaction(s) of a type this library does not model"},
      {diag.unresolvedYieldDelegations, "delegated fission-yield table(s) that did not resolve"},
      {diag.unparsedReactionTargets, "reaction target(s) that did not parse (the product lost)"},
  };
  std::fprintf(stderr, "cram: WARNING - %s: unusable entries dropped:\n", what.c_str());
  for (const auto& [n, text] : counts) {
    if (n > 0)
      std::fprintf(stderr, "cram:   %d %s\n", n, text);
  }
}

std::vector<ChainReaction> loadDocument(DepletionChain& chain, const pugi::xml_document& doc,
                                        const std::string& what, ChainXmlDiagnostics* diagnostics) {
  pugi::xml_node root = doc.child("depletion_chain");
  if (!root)
    throw std::runtime_error("cram: " + what + " has no <depletion_chain> root element");

  ChainXmlDiagnostics diag;

  // First pass: register every nuclide so decay/reaction targets resolve, and
  // index the elements by name so a yield delegation can find its source
  // whether it appears before or after the nuclide that delegates to it.
  std::unordered_map<std::string, pugi::xml_node> byName;
  for (pugi::xml_node nuc : root.children("nuclide")) {
    const char* name = nuc.attribute("name").as_string();
    if (const std::optional<Zai> z = parseNuclideName(name)) {
      chain.add(*z);
      byName.emplace(name, nuc);
    } else {
      ++diag.unparsedNuclides;
    }
  }

  // Second pass: decay, fission yields, reaction topology.
  std::vector<ChainReaction> reactions;
  std::vector<std::pair<Zai, std::string>> delegated;
  for (pugi::xml_node nuc : root.children("nuclide")) {
    const std::optional<Zai> z = parseNuclideName(nuc.attribute("name").as_string());
    if (!z)
      continue;  // already counted
    parseDecay(nuc, *z, chain, diag);
    if (pugi::xml_node nfy = nuc.child("neutron_fission_yields")) {
      const std::string_view from = nfy.attribute("parent").as_string();
      if (from.empty())
        parseFissionYields(nfy, *z, chain, diag);
      else
        delegated.emplace_back(*z, std::string(from));  // resolved below
    }
    parseReactions(nuc, *z, chain, reactions, diag);
  }

  // Third pass: the delegating nuclides. Deferred to here because the nuclide
  // whose tables they borrow may itself delegate, and may appear anywhere in
  // the file.
  for (const auto& [z, from] : delegated) {
    pugi::xml_node nfy = resolveYieldDelegation(byName, from);
    if (nfy)
      parseFissionYields(nfy, z, chain, diag);
    else
      ++diag.unresolvedYieldDelegations;
  }

  // Reported through the out-parameter when the caller asked for it, and only
  // written to stderr otherwise -- the same convention as decayMatrix().
  if (diagnostics != nullptr)
    *diagnostics = diag;
  else
    warnDiagnostics(diag, what);
  return reactions;
}

}  // namespace

std::vector<ChainReaction> loadDepletionChainXml(DepletionChain& chain, const std::string& path,
                                                 ChainXmlDiagnostics* diagnostics) {
  pugi::xml_document doc;
  pugi::xml_parse_result res = doc.load_file(path.c_str());
  if (!res)
    throw std::runtime_error("cram: failed to parse depletion_chain XML '" + path +
                             "': " + res.description());
  return loadDocument(chain, doc, "'" + path + "'", diagnostics);
}

std::vector<ChainReaction> loadDepletionChainXmlString(DepletionChain& chain, std::string_view xml,
                                                       ChainXmlDiagnostics* diagnostics) {
  pugi::xml_document doc;
  pugi::xml_parse_result res = doc.load_buffer(xml.data(), xml.size());
  if (!res)
    throw std::runtime_error(std::string("cram: failed to parse depletion_chain XML: ") +
                             res.description());
  return loadDocument(chain, doc, "in-memory depletion_chain XML", diagnostics);
}

#else  // !CRAM_WITH_CHAIN_XML

std::vector<ChainReaction> loadDepletionChainXml(DepletionChain&, const std::string& path,
                                                 ChainXmlDiagnostics*) {
  throw std::runtime_error("cram: cannot read '" + path +
                           "': built without CRAM_WITH_CHAIN_XML (the OpenMC chain reader)");
}

std::vector<ChainReaction> loadDepletionChainXmlString(DepletionChain&, std::string_view,
                                                       ChainXmlDiagnostics*) {
  throw std::runtime_error("cram: built without CRAM_WITH_CHAIN_XML (the OpenMC chain reader)");
}

#endif

}  // namespace cram
