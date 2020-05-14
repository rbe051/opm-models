// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup multiDomain
 * \brief A helper classe that converts a Dune::MultiTypeBlockMatrix into a plain Dune::BCRSMatrix
 * This code is mostly borrowed from the DuMuX multidomain module: https://dumux.org/
 */
#ifndef OPM_MATRIX_CONVERTER
#define OPM_MATRIX_CONVERTER

#include <cmath>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <chrono>

namespace Opm
{

/*!
 * \ingroup multiDomain
 * \brief A helper classe that converts a Dune::MultiTypeBlockMatrix into a plain Dune::BCRSMatrix
 * TODO: allow block sizes for BCRSMatrix other than 1x1 ?
 *
 */
template <class MultiTypeBlockMatrix, class MatrixBlock, class Scalar = double>
class MatrixConverter
{
    using BCRSMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

public:
    /*!
     * \brief Converts the matrix to a type the IterativeSolverBackend can handle
     *
     * \param A The original multitype blockmatrix
     */
    static auto multiTypeToBCRSMatrix(const MultiTypeBlockMatrix &A)
    {
        // get the size for the converted matrix
        const auto numRows = getNumRows_(A);

        // create an empty BCRS matrix with 1x1 blocks
        auto M = BCRSMatrix(numRows, numRows, BCRSMatrix::random);

        // set the occupation pattern and copy the values
        setOccupationPattern_(M, A);
        copyValues_(M, A);

        return M;
    }

private:
    /*!
     * \brief Sets the occupation pattern and indices for the converted matrix
     *
     * \param M The converted matrix
     * \param A The original multitype blockmatrix
     */
    static void setOccupationPattern_(BCRSMatrix &M, const MultiTypeBlockMatrix &A)
    {
        // prepare the occupation pattern
        const auto numRows = M.N();
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numRows, numRows);

        // lambda function to fill the occupation pattern
        auto addIndices = [&occupationPattern](const auto &subMatrix, const std::size_t startRow, const std::size_t startCol) {
            using std::abs;
            static const Scalar eps = -1.0;

            using BlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
            const auto blockSizeI = BlockType::rows;
            const auto blockSizeJ = BlockType::cols;

            const auto copyToScalar = MatrixBlock::rows == 1 && MatrixBlock::cols == 1;
            if (!copyToScalar && (blockSizeI != MatrixBlock::rows || blockSizeJ != MatrixBlock::rows))
                throw std::logic_error("Matrix block size must match to merge them");

            for (auto row = subMatrix.begin(); row != subMatrix.end(); ++row)
            {
                for (auto col = row->begin(); col != row->end(); ++col)
                {
                    if (copyToScalar)
                    {
                        for (std::size_t i = 0; i < blockSizeI; ++i)
                            for (std::size_t j = 0; j < blockSizeJ; ++j)
                                if (abs(subMatrix[row.index()][col.index()][i][j]) > eps)
                                    occupationPattern.add(startRow + row.index() * blockSizeI + i, startCol + col.index() * blockSizeJ + j);
                    }
                    else if (subMatrix[row.index()][col.index()].infinity_norm() > eps)
                        occupationPattern.add(startRow + row.index(), startCol + col.index());
                }
            }
        };

        // fill the pattern
        std::size_t rowIndex = 0;
        const auto copyToScalar = MatrixBlock::rows == 1 && MatrixBlock::cols == 1;

        Dune::Hybrid::forEach(A, [&addIndices, &rowIndex, numRows, copyToScalar](const auto &rowOfMultiTypeMatrix) {
            std::size_t colIndex = 0;
            Dune::Hybrid::forEach(rowOfMultiTypeMatrix, [&addIndices, &colIndex, &rowIndex, numRows, copyToScalar](const auto &subMatrix) {
                addIndices(subMatrix, rowIndex, colIndex);

                using SubBlockType = typename std::decay_t<decltype(subMatrix)>::block_type;

                if (copyToScalar)
                    colIndex += SubBlockType::cols * subMatrix.M();
                else
                    colIndex += subMatrix.M();

                // if we have arrived at the right side of the matrix, increase the row index
                if (colIndex == numRows)
                {
                    if (copyToScalar)
                        rowIndex += SubBlockType::rows * subMatrix.N();
                    else
                        rowIndex += subMatrix.N();
                }
            });
        });

        occupationPattern.exportIdx(M);
    }

    /*!
     * \brief Sets the occupation pattern (i.e. indices) for the converted matrix
     *
     * \param M The converted matrix
     * \param A The original subMatrix of the multitype blockmatrix
     */
    static void
    copyValues_(BCRSMatrix &M, const MultiTypeBlockMatrix &A)
    {
        // get number of rows
        const auto numRows = M.N();
        const auto numCols = M.M();

        // lambda function to copy the values
        auto copyValues = [&M](const auto &subMatrix, const std::size_t startRow, const std::size_t startCol) {
            using std::abs;
            static const Scalar eps = -1.0;

            using BlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
            const auto blockSizeI = BlockType::rows;
            const auto blockSizeJ = BlockType::cols;

            const auto copyToScalar = MatrixBlock::rows == 1 && MatrixBlock::cols == 1;
            if (!copyToScalar && (blockSizeI != MatrixBlock::rows || blockSizeJ != MatrixBlock::rows))
                throw std::logic_error("Matrix block size must match to merge them");

            for (auto row = subMatrix.begin(); row != subMatrix.end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    if (copyToScalar)
                    {
                        for (std::size_t i = 0; i < blockSizeI; ++i)
                            for (std::size_t j = 0; j < blockSizeJ; ++j)
                                if (abs(subMatrix[row.index()][col.index()][i][j]) > eps)
                                    M[startRow + row.index() * blockSizeI + i][startCol + col.index() * blockSizeJ + j] = subMatrix[row.index()][col.index()][i][j];
                    }
                    else
                    {
                        if (subMatrix[row.index()][col.index()].infinity_norm() > eps)
                            M[startRow + row.index()][startCol + col.index()] = subMatrix[row.index()][col.index()];
                    }
        };

        std::size_t rowIndex = 0;

        const auto copyToScalar = MatrixBlock::rows == 1 && MatrixBlock::cols == 1;

        Dune::Hybrid::forEach(A, [&copyValues, &rowIndex, numRows, copyToScalar](const auto &rowOfMultiTypeMatrix) {
            std::size_t colIndex = 0;
            Dune::Hybrid::forEach(rowOfMultiTypeMatrix, [&copyValues, &colIndex, &rowIndex, numRows, copyToScalar](const auto &subMatrix) {
                copyValues(subMatrix, rowIndex, colIndex);

                using SubBlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
                if (copyToScalar)
                    colIndex += SubBlockType::cols * subMatrix.M();
                else
                    colIndex += subMatrix.M();

                // if we have arrived at the right side of the matrix, increase the row index
                if (colIndex == numRows)
                {
                    if (copyToScalar)
                        rowIndex += SubBlockType::rows * subMatrix.N();
                    else
                        rowIndex += subMatrix.N();
                }
            });
        });
    }

    /*!
     * \brief Calculates the total number of rows (== number of cols) for the converted matrix
     *
     * \param A The original multitype blockmatrix
     */
    static std::size_t getNumRows_(const MultiTypeBlockMatrix &A)
    {
        // iterate over the first row of the multitype blockmatrix
        std::size_t numRows = 0;
        const auto copyToScalar = MatrixBlock::rows == 1 && MatrixBlock::cols == 1;

        Dune::Hybrid::forEach(Dune::Hybrid::elementAt(A, Dune::Indices::_0), [&numRows, copyToScalar](const auto &subMatrixInFirstRow) {
            // the number of cols of the individual submatrice's block equals the respective number of equations.
            if (copyToScalar)
            {
                const auto numEq = std::decay_t<decltype(subMatrixInFirstRow)>::block_type::cols;
                numRows += numEq * subMatrixInFirstRow.M();
            }
            else
                numRows += subMatrixInFirstRow.M();
        });

        return numRows;
    }
};

/*!
 * \ingroup multiDomain
 * \brief A helper classe that converts a Dune::MultiTypeBlockVector into a plain Dune::BlockVector and transfers back values
 */
template <class MultiTypeBlockVector, class VectorBlock, class Scalar = double>
class VectorConverter
{
    using BlockVector = typename Dune::BlockVector<VectorBlock>;

public:
    /*!
     * \brief Converts a Dune::MultiTypeBlockVector to a plain 1x1 Dune::BlockVector
     *
     * \param b The original multitype blockvector
     */
    static auto multiTypeToBlockVector(const MultiTypeBlockVector &b)
    {
        const auto size = getSize_(b);

        BlockVector bTmp;
        bTmp.resize(size);

        const auto copyToScalar = VectorBlock::size() == 1;

        std::size_t startIndex = 0;
        Dune::Hybrid::forEach(b, [&bTmp, &startIndex, copyToScalar](const auto &subVector) {
            const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();

            if (!copyToScalar && numEq != VectorBlock::size())
                throw std::logic_error("Vector block sizes must match to merge them");

            for (std::size_t i = 0; i < subVector.size(); ++i)
                if (copyToScalar)
                {
                    for (std::size_t j = 0; j < numEq; ++j)
                        bTmp[startIndex + i * numEq + j] = subVector[i][j];
                }
                else
                    bTmp[startIndex + i] = subVector[i];

            if (copyToScalar)
                startIndex += numEq * subVector.size();
            else
                startIndex += subVector.size();
        });

        return bTmp;
    }

    /*!
     * \brief Copys the entries of a Dune::BlockVector to a Dune::MultiTypeBlockVector
     *
     * \param x The multitype blockvector where the values are copied to
     * \param y The regular blockvector where the values are copied from
     */
    static void retrieveValues(MultiTypeBlockVector &x, const BlockVector &y)
    {
        std::size_t startIndex = 0;

        Dune::Hybrid::forEach(x, [&y, &startIndex](auto &subVector) {
            const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();
            const auto copyToScalar = VectorBlock::size() == 1;
            if (!copyToScalar && numEq != VectorBlock::size())
                throw std::logic_error("Vector block sizes must match to merge them");

            for (std::size_t i = 0; i < subVector.size(); ++i)
                if (copyToScalar)
                {
                    for (std::size_t j = 0; j < numEq; ++j){
                        throw std::logic_error("Not implemented");
                        //subVector[i][j] = y[startIndex + i * numEq + j];
                    }
                }
                else
                    subVector[i] = y[startIndex + i];

            if (copyToScalar)
                startIndex += numEq * subVector.size();
            else
                startIndex += subVector.size();
        });
    }

private:
    /*!
     * \brief Returns the size of the expanded multitype block vector
     *
     * \param b The multitype blockvector
     */
    static std::size_t getSize_(const MultiTypeBlockVector &b)
    {
        std::size_t size = 0;
        const auto copyToScalar = VectorBlock::size() == 1;

        Dune::Hybrid::forEach(b, [&size, copyToScalar](const auto &subVector) {
            // the size of the individual vector blocks equals the respective number of equations.
            if (copyToScalar)
            {
                const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();
                size += numEq * subVector.size();
            }
            else
                size += subVector.size();
        });
        return size;
    }
};

} // end namespace Opm

#endif

