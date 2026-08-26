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
    if (type == "sf") {
      // Products from the fission-yield table, if present. rtyp 6.0 is the
      // ENDF code for spontaneous fission; decayDaughter() never consults it
      // for a fission mode, but it keeps the mode self-describing.
      d.modes.push_back(
          DecayMode{.rtyp = 6.0, .branching = branching, .finalState = 0, .isFission = true});
      continue;
    }
    const std::optional<Zai> daughter = parseNuclideName(dn.attribute("target").as_string());
    if (!daughter) {
      ++diag.unparsedDecayTargets;  // mode dropped; its branching is lost
      continue;
    }
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

void parseFissionYields(const pugi::xml_node& nuc, const Zai& parent, DepletionChain& chain,
                        ChainXmlDiagnostics& diag) {
  pugi::xml_node nfy = nuc.child("neutron_fission_yields");
  if (!nfy)
    return;
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

void parseReactions(const pugi::xml_node& nuc, const Zai& parent, std::vector<ChainReaction>& out,
                    ChainXmlDiagnostics& diag) {
  for (pugi::xml_node rn : nuc.children("reaction")) {
    const std::optional<ReactionType> type = reactionTypeFromName(rn.attribute("type").as_string());
    if (!type) {
      ++diag.unmodeledReactions;
      continue;
    }
    const std::string_view target = rn.attribute("target").as_string();
    std::optional<Zai> product;
    if (!target.empty() && target != "Nothing")
      product = parseNuclideName(target);
    out.push_back(ChainReaction{.parent = parent,
                                .type = *type,
                                .target = product,
                                .q = rn.attribute("Q").as_double(0.0),
                                .branching = rn.attribute("branching_ratio").as_double(1.0)});
  }
}

std::vector<ChainReaction> loadDocument(DepletionChain& chain, const pugi::xml_document& doc,
                                        const std::string& what, ChainXmlDiagnostics* diagnostics) {
  pugi::xml_node root = doc.child("depletion_chain");
  if (!root)
    throw std::runtime_error("cram: " + what + " has no <depletion_chain> root element");

  ChainXmlDiagnostics diag;

  // First pass: register every nuclide so decay/reaction targets resolve.
  for (pugi::xml_node nuc : root.children("nuclide")) {
    if (const std::optional<Zai> z = parseNuclideName(nuc.attribute("name").as_string()))
      chain.add(*z);
    else
      ++diag.unparsedNuclides;
  }

  // Second pass: decay, fission yields, reaction topology.
  std::vector<ChainReaction> reactions;
  for (pugi::xml_node nuc : root.children("nuclide")) {
    const std::optional<Zai> z = parseNuclideName(nuc.attribute("name").as_string());
    if (!z)
      continue;  // already counted
    parseDecay(nuc, *z, chain, diag);
    parseFissionYields(nuc, *z, chain, diag);
    parseReactions(nuc, *z, reactions, diag);
  }

  // Reported through the out-parameter when the caller asked for it, and only
  // written to stderr otherwise -- the same convention as decayMatrix().
  if (diagnostics != nullptr) {
    *diagnostics = diag;
  } else if (!diag.clean()) {
    std::fprintf(stderr,
                 "cram: WARNING - %s: dropped %d nuclide(s), %d decay target(s), %d yield "
                 "product(s) that did not parse and %d reaction(s) of unmodeled type.\n",
                 what.c_str(), diag.unparsedNuclides, diag.unparsedDecayTargets,
                 diag.unparsedYieldProducts, diag.unmodeledReactions);
  }
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
