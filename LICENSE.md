# License and attributions

## This project

Copyright (c) 2026. Released under the BSD 3-Clause License (below). You may
adapt this header with your name/organization as appropriate.

## Third-party material bundled or derived here

The following components carry their own licenses and attribution obligations.
They are noted because they were copied or derived from external sources.

1. CRAM coefficients (`src/cram.cpp`)
   The IPF-CRAM-16 and CRAM-48 pole/residue tables were transcribed from
   OpenMC's `openmc/deplete/cram.py`, which is distributed under the MIT
   License (Copyright (c) Massachusetts Institute of Technology, UChicago
   Argonne LLC, and OpenMC contributors). The underlying coefficients are from
   M. Pusa, "Higher-Order Chebyshev Rational Approximation Method and
   Application to Burnup Equations," Nucl. Sci. Eng. 182:3, 297-318 (2016).

2. Sample ENDF data (`tests/integration/data/`)
   * `nfy-Pu239.endf` is copied from the njoy/ENDFtk repository's example
     resources (BSD 3-Clause, Copyright (c) njoy contributors).
   * `decay_am242m.endf` was assembled from an MF8/MT457 record used in
     ENDFtk's unit-test fixtures (same source/license). These files are
     evaluated nuclear data included solely for testing the reader.

3. ENDFtk (optional dependency, fetched at build time with -DWITH_ENDFTK=ON)
   njoy/ENDFtk and its dependencies (njoy/tools, spdlog, fast_float) are not
   bundled here; they are fetched via CMake FetchContent and retain their own
   licenses (BSD 3-Clause / MIT).

-------------------------------------------------------------------------------
BSD 3-Clause License

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
