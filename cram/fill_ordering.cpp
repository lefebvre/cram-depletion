#include <Eigen/OrderingMethods>
#include <algorithm>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "cram/fill_ordering_internal.hpp"

namespace cram::internal {
namespace {

// Strongly connected components of the graph with an edge col -> row for every
// stored off-diagonal entry, in topological order of the condensation: a
// component comes before every component it has an edge into. Tarjan's
// algorithm, iterative so a long chain cannot exhaust the stack; it emits
// components sinks-first, so the result is reversed at the end.
template <typename Scalar>
std::vector<std::vector<int>> componentsInTopologicalOrder(const Eigen::SparseMatrix<Scalar>& A) {
  const int n = static_cast<int>(A.rows());
  std::vector<std::vector<int>> out(static_cast<std::size_t>(n));
  for (int c = 0; c < A.outerSize(); ++c)
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, c); it; ++it)
      if (it.row() != it.col())
        out[static_cast<std::size_t>(it.col())].push_back(static_cast<int>(it.row()));

  std::vector<int> index(static_cast<std::size_t>(n), -1);
  std::vector<int> low(static_cast<std::size_t>(n), 0);
  std::vector<char> onStack(static_cast<std::size_t>(n), 0);
  std::vector<int> stack;
  std::vector<std::vector<int>> components;
  // (vertex, next edge to visit) frames standing in for the recursion.
  std::vector<std::pair<int, std::size_t>> frames;
  int counter = 0;
  for (int root = 0; root < n; ++root) {
    if (index[static_cast<std::size_t>(root)] >= 0)
      continue;
    frames.emplace_back(root, 0);
    while (!frames.empty()) {
      auto& [v, edge] = frames.back();
      const auto vi = static_cast<std::size_t>(v);
      if (edge == 0 && index[vi] < 0) {
        index[vi] = low[vi] = counter++;
        stack.push_back(v);
        onStack[vi] = 1;
      }
      if (edge < out[vi].size()) {
        const int w = out[vi][edge++];
        const auto wi = static_cast<std::size_t>(w);
        if (index[wi] < 0)
          frames.emplace_back(w, 0);  // invalidates v/edge; re-read next pass
        else if (onStack[wi] != 0)
          low[vi] = std::min(low[vi], index[wi]);
        continue;
      }
      if (low[vi] == index[vi]) {
        std::vector<int> comp;
        int w = -1;
        do {
          w = stack.back();
          stack.pop_back();
          onStack[static_cast<std::size_t>(w)] = 0;
          comp.push_back(w);
        } while (w != v);
        components.push_back(std::move(comp));
      }
      const int finished = v;
      frames.pop_back();
      if (!frames.empty()) {
        const auto pi = static_cast<std::size_t>(frames.back().first);
        low[pi] = std::min(low[pi], low[static_cast<std::size_t>(finished)]);
      }
    }
  }
  std::reverse(components.begin(), components.end());
  return components;
}

}  // namespace

template <typename Scalar>
SymmetricPermutation fillReducingOrder(const Eigen::SparseMatrix<Scalar>& A) {
  const int n = static_cast<int>(A.rows());
  SymmetricPermutation perm(n);
  std::vector<int> local(static_cast<std::size_t>(n), -1);
  int next = 0;
  for (std::vector<int>& comp : componentsInTopologicalOrder(A)) {
    // A block of one or two nuclides has nothing for COLAMD to improve.
    if (comp.size() < 3) {
      for (int v : comp)
        perm.indices()(v) = next++;
      continue;
    }
    // COLAMD on the block's own pattern, applied to rows and columns alike so
    // the block's diagonal stays on the diagonal, where the pivots are.
    std::sort(comp.begin(), comp.end());
    const auto m = static_cast<int>(comp.size());
    for (int q = 0; q < m; ++q)
      local[static_cast<std::size_t>(comp[static_cast<std::size_t>(q)])] = q;
    std::vector<Eigen::Triplet<double>> t;
    for (int q = 0; q < m; ++q) {
      t.emplace_back(q, q, 1.0);
      for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(
               A, comp[static_cast<std::size_t>(q)]);
           it; ++it) {
        const int r = local[static_cast<std::size_t>(it.row())];
        if (r >= 0)
          t.emplace_back(r, q, 1.0);
      }
    }
    Eigen::SparseMatrix<double> block(m, m);
    block.setFromTriplets(t.begin(), t.end());
    block.makeCompressed();
    SymmetricPermutation blockPerm;
    Eigen::COLAMDOrdering<int>()(block, blockPerm);
    for (int q = 0; q < m; ++q) {
      perm.indices()(comp[static_cast<std::size_t>(q)]) = next + blockPerm.indices()(q);
      local[static_cast<std::size_t>(comp[static_cast<std::size_t>(q)])] = -1;
    }
    next += m;
  }
  return perm;
}

template SymmetricPermutation fillReducingOrder(const Eigen::SparseMatrix<double>&);
template SymmetricPermutation fillReducingOrder(const Eigen::SparseMatrix<std::complex<double>>&);

}  // namespace cram::internal
