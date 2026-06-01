#include "cram/chain_xml.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

#include <pugixml.hpp>

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
    if (kSymbols[z] == sym)
      return z;
  return -1;
}

ReactionType reactionTypeFromString(std::string_view s, bool& known) {
  known = true;
  if (s == "fission")
    return ReactionType::Fission;
  if (s == "(n,gamma)")
    return ReactionType::NGamma;
  if (s == "(n,2n)")
    return ReactionType::N2n;
  if (s == "(n,3n)")
    return ReactionType::N3n;
  if (s == "(n,4n)")
    return ReactionType::N4n;
  if (s == "(n,a)")
    return ReactionType::NAlpha;
  if (s == "(n,p)")
    return ReactionType::NProton;
  known = false;
  return ReactionType::NGamma;
}

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

void parseDecay(const pugi::xml_node& nuc, const Zai& parent, DepletionChain& chain) {
  const double halfLife = nuc.attribute("half_life").as_double(0.0);
  if (halfLife <= 0.0)
    return;  // stable (or no half-life given)

  DecayData d;
  d.halfLife = halfLife;
  for (pugi::xml_node dn : nuc.children("decay")) {
    DecayMode m;
    m.branching = dn.attribute("branching_ratio").as_double(0.0);
    const std::string type = dn.attribute("type").as_string();
    if (type == "sf") {
      m.isFission = true;  // products from the fission-yield table, if present
    } else {
      m.hasDaughter = true;  // {0,0,0} (dropped) if the target does not parse
      parseNuclideName(dn.attribute("target").as_string(), m.daughter);
    }
    d.modes.push_back(m);
  }
  chain.setDecay(parent, std::move(d));
}

void parseFissionYields(const pugi::xml_node& nuc, const Zai& parent, DepletionChain& chain) {
  pugi::xml_node nfy = nuc.child("neutron_fission_yields");
  if (!nfy)
    return;
  for (pugi::xml_node fy : nfy.children("fission_yields")) {
    FissionYields y;
    y.energy = fy.attribute("energy").as_double(0.0);
    const auto products = tokenize(fy.child("products").text().get());
    const auto data = tokenize(fy.child("data").text().get());
    const std::size_t n = std::min(products.size(), data.size());
    for (std::size_t k = 0; k < n; ++k) {
      Zai prod;
      if (!parseNuclideName(products[k], prod))
        continue;
      y.products.emplace_back(prod, std::stod(data[k]));
    }
    chain.addFissionYields(parent, std::move(y));
  }
}

void parseReactions(const pugi::xml_node& nuc, const Zai& parent, std::vector<ChainReaction>& out) {
  for (pugi::xml_node rn : nuc.children("reaction")) {
    bool known = false;
    ReactionType type = reactionTypeFromString(rn.attribute("type").as_string(), known);
    if (!known)
      continue;
    ChainReaction r;
    r.parent = parent;
    r.type = type;
    r.q = rn.attribute("Q").as_double(0.0);
    r.branching = rn.attribute("branching_ratio").as_double(1.0);
    const std::string target = rn.attribute("target").as_string();
    if (!target.empty() && target != "Nothing")
      r.hasTarget = parseNuclideName(target, r.target);
    out.push_back(r);
  }
}

}  // namespace

std::string elementSymbol(int z) {
  if (z < 1 || z >= static_cast<int>(kSymbols.size()))
    return {};
  return std::string(kSymbols[z]);
}

bool parseNuclideName(const std::string& name, Zai& out) {
  std::size_t i = 0;
  while (i < name.size() && (std::isalpha(static_cast<unsigned char>(name[i])) != 0))
    ++i;
  if (i == 0 || i >= name.size())
    return false;
  const int z = symbolToZ(std::string_view(name).substr(0, i));
  if (z < 0)
    return false;

  // Mass number digits follow the symbol.
  std::size_t j = i;
  while (j < name.size() && (std::isdigit(static_cast<unsigned char>(name[j])) != 0))
    ++j;
  if (j == i)
    return false;
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

  out = Zai{z, a, iso};
  return true;
}

std::vector<ChainReaction> loadDepletionChainXml(DepletionChain& chain, const std::string& path) {
  pugi::xml_document doc;
  pugi::xml_parse_result res = doc.load_file(path.c_str());
  if (!res)
    throw std::runtime_error("cram: failed to parse depletion_chain XML '" + path +
                             "': " + res.description());

  pugi::xml_node root = doc.child("depletion_chain");
  if (!root)
    throw std::runtime_error("cram: '" + path + "' has no <depletion_chain> root element");

  // First pass: register every nuclide so decay/reaction targets resolve.
  std::vector<ChainReaction> reactions;
  for (pugi::xml_node nuc : root.children("nuclide")) {
    Zai z;
    if (parseNuclideName(nuc.attribute("name").as_string(), z))
      chain.add(z);
  }

  // Second pass: decay, fission yields, reaction topology.
  for (pugi::xml_node nuc : root.children("nuclide")) {
    Zai z;
    if (!parseNuclideName(nuc.attribute("name").as_string(), z))
      continue;
    parseDecay(nuc, z, chain);
    parseFissionYields(nuc, z, chain);
    parseReactions(nuc, z, reactions);
  }
  return reactions;
}

}  // namespace cram
