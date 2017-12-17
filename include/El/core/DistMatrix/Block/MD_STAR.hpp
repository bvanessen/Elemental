/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DISTMATRIX_BLOCKCYCLIC_MD_STAR_HPP
#define EL_DISTMATRIX_BLOCKCYCLIC_MD_STAR_HPP

#include "El/core/DistMatrix_decl.hpp"
#include "El/core/Grid.hpp"
#include "El/macros.h"

namespace El
{

// Partial specialization to A[MD,* ].
//
// The columns of these distributed matrices will be distributed like
// "Matrix Diagonals" (MD). It is important to recognize that the diagonal
// of a sufficiently large distributed matrix is distributed amongst the
// entire process grid if and only if the dimensions of the process grid
// are coprime.
template<typename Ring>
class DistMatrix<Ring,Dist::MD,Dist::STAR,DistWrap::BLOCK> : public BlockMatrix<Ring>
{
public:
    typedef AbstractDistMatrix<Ring> absType;
    typedef BlockMatrix<Ring> blockType;
    typedef DistMatrix<Ring,Dist::MD,Dist::STAR,DistWrap::BLOCK> type;
    typedef DistMatrix<Ring,Dist::STAR,Dist::MD,DistWrap::BLOCK> transType;
    typedef DistMatrix<Ring,Dist::MD,Dist::STAR,DistWrap::BLOCK> diagType;

    // Constructors and destructors
    // ============================

    // Create a 0 x 0 distributed matrix with default (and unpinned) block size
    DistMatrix(const El::Grid& grid=Grid::Default(), int root=0);

    // Create a 0 x 0 distributed matrix with fixed block size
    DistMatrix
    (const El::Grid& grid, Int blockHeight, Int blockWidth, int root=0);

    // Create a height x width distributed matrix with default block size
    DistMatrix
    (Int height, Int width, const El::Grid& grid=Grid::Default(), int root=0);

    // Create a height x width distributed matrix with fixed block size
    DistMatrix
    (Int height, Int width, const El::Grid& grid,
      Int blockHeight, Int blockWidth, int root=0 );

    // Create a copy of distributed matrix A (redistributing if necessary)
    DistMatrix(const type& A);
    DistMatrix(const absType& A );
    DistMatrix(const blockType& A);
    template<Dist colDist,Dist rowDist>
    DistMatrix(const DistMatrix<Ring,colDist,rowDist,DistWrap::BLOCK>& A );
    template<Dist colDist,Dist rowDist>
    DistMatrix(const DistMatrix<Ring,colDist,rowDist,DistWrap::ELEMENT>& A);

    // Move constructor
    DistMatrix(type&& A) EL_NO_EXCEPT;

    // Destructor
    ~DistMatrix();

    type* Copy() const override;
    type* Construct
    (const El::Grid& grid, int root) const override;
    transType* ConstructTranspose
    (const El::Grid& grid, int root) const override;
    diagType* ConstructDiagonal
    (const El::Grid& grid, int root) const override;

    // Operator overloading
    // ====================

    // Return a view of a contiguous submatrix
    // ---------------------------------------
          type operator()(Range<Int> I, Range<Int> J );
    const type operator()(Range<Int> I, Range<Int> J ) const;

    // Return a copy of a (generally non-contiguous) submatrix
    // -------------------------------------------------------
    type operator()(Range<Int> I, const std::vector<Int>& J) const;
    type operator()(const std::vector<Int>& I, Range<Int> J ) const;
    type operator()(const std::vector<Int>& I, const std::vector<Int>& J ) const;

    // Make a copy
    // -----------
    type& operator=(const absType& A );
    type& operator=(const blockType& A );
    type& operator=(const DistMatrix<Ring,Dist::MC,  Dist::MR  ,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::MC,  Dist::STAR,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::MR  ,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::MD,  Dist::STAR,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::MD  ,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::MR,  Dist::MC  ,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::MR,  Dist::STAR,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::MC  ,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::VC,  Dist::STAR,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::VC  ,DistWrap::BLOCK>& A);
    type& operator=(const DistMatrix<Ring,Dist::VR,  Dist::STAR,DistWrap::BLOCK>& A );
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::VR  ,DistWrap::BLOCK>& A);
    type& operator=(const DistMatrix<Ring,Dist::STAR,Dist::STAR,DistWrap::BLOCK>& A);
    type& operator=(const DistMatrix<Ring,Dist::CIRC,Dist::CIRC,DistWrap::BLOCK>& A);
    template<Dist colDist,Dist rowDist>
    type& operator=(const DistMatrix<Ring,colDist,rowDist,DistWrap::ELEMENT>& A);

    // Move assignment
    // ---------------
    type& operator=(type&& A);

    // Rescaling
    // ---------
    const type& operator*=(Ring alpha );

    // Addition/subtraction
    // --------------------
    const type& operator+=(const blockType& A);
    const type& operator+=(const absType& A);
    const type& operator-=(const blockType& A);
    const type& operator-=(const absType& A);

    // Basic queries
    // =============
    Dist ColDist()             const override EL_NO_EXCEPT;
    Dist RowDist()             const override EL_NO_EXCEPT;
    Dist PartialColDist()      const override EL_NO_EXCEPT;
    Dist PartialRowDist()      const override EL_NO_EXCEPT;
    Dist PartialUnionColDist() const override EL_NO_EXCEPT;
    Dist PartialUnionRowDist() const override EL_NO_EXCEPT;
    Dist CollectedColDist()    const override EL_NO_EXCEPT;
    Dist CollectedRowDist()    const override EL_NO_EXCEPT;

    mpi::Comm DistComm()            const EL_NO_EXCEPT override;
    mpi::Comm CrossComm()           const EL_NO_EXCEPT override;
    mpi::Comm RedundantComm()       const EL_NO_EXCEPT override;
    mpi::Comm ColComm()             const EL_NO_EXCEPT override;
    mpi::Comm RowComm()             const EL_NO_EXCEPT override;
    mpi::Comm PartialColComm()      const EL_NO_EXCEPT override;
    mpi::Comm PartialRowComm()      const EL_NO_EXCEPT override;
    mpi::Comm PartialUnionColComm() const EL_NO_EXCEPT override;
    mpi::Comm PartialUnionRowComm() const EL_NO_EXCEPT override;

    int DistSize()              const EL_NO_EXCEPT override;
    int CrossSize()             const EL_NO_EXCEPT override;
    int RedundantSize()         const EL_NO_EXCEPT override;
    int ColStride()             const EL_NO_EXCEPT override;
    int RowStride()             const EL_NO_EXCEPT override;
    int PartialColStride()      const EL_NO_EXCEPT override;
    int PartialRowStride()      const EL_NO_EXCEPT override;
    int PartialUnionColStride() const EL_NO_EXCEPT override;
    int PartialUnionRowStride() const EL_NO_EXCEPT override;

    int DistRank()            const EL_NO_EXCEPT override;
    int CrossRank()           const EL_NO_EXCEPT override;
    int RedundantRank()       const EL_NO_EXCEPT override;
    int ColRank()             const EL_NO_EXCEPT override;
    int RowRank()             const EL_NO_EXCEPT override;
    int PartialColRank()      const EL_NO_EXCEPT override;
    int PartialRowRank()      const EL_NO_EXCEPT override;
    int PartialUnionColRank() const EL_NO_EXCEPT override;
    int PartialUnionRowRank() const EL_NO_EXCEPT override;

    template<typename S,Dist U,Dist V,DistWrap wrap> friend class DistMatrix;
};

} // namespace El

#endif // ifndef EL_DISTMATRIX_BLOCKCYCLIC_MD_STAR_HPP
