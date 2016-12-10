//---------------------------------------------------------------------------
//    $Id: solver_utils.cc 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2008 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <solvers/solver_utils.h>

//---------------------------------------------------------------------------
void
SolverUtils::check_diagonal(const BlockSparseMatrix<double>& A)
{
  // Find zero element on the diagonal
  std::vector<unsigned int> zero_diag;
  for (unsigned int i=0; i < A.m(); i++)
    {
      const double diag_element = A.diag_element(i);
      if (fabs(diag_element) == 0.0)
	zero_diag.push_back(i);
    }
  // Give a warning:
  Assert(zero_diag.size() != 0,
	 ExcMessage("WARNING: Zeros in the diagonal"));
  
  if (zero_diag.size() > 0)
    FcstUtilities::log<<"WARNING: "<<zero_diag.size()<<" diagonal elements are zero"<<std::endl;
}
//---------------------------------------------------------------------------
void
SolverUtils::output_diagonal(const BlockSparseMatrix<double>& A)
{
  std::ofstream file;
  FcstUtilities::log<<"============================"<<std::endl;
  FcstUtilities::log<<"== DIAGONAL TERMS =="<<std::endl;
  for (unsigned int i = 0; i< A.m(); ++i)
    FcstUtilities::log<<"Diagonal element in row "<<i<<" is "<<A.diag_element(i)<<std::endl;
  FcstUtilities::log<<"============================"<<std::endl;
}
//---------------------------------------------------------------------------
void
SolverUtils::print_diagonal(const BlockSparseMatrix<double>& A, 
                            const std::string& in_file)
{
  std::ofstream file;
  file.open(in_file.c_str());//"diag_matrix.dat");
  file<<"============================"<<std::endl;
  file<<"== DIAGONAL TERMS =="<<std::endl;
  for (unsigned int i = 0; i< A.m(); ++i)
    file<<"Diagonal element in row "<<i<<" is "<<A.diag_element(i)<<std::endl;
  file<<"============================"<<std::endl;
  file.close();
}

//---------------------------------------------------------------------------
void 
SolverUtils::repair_diagonal(BlockSparseMatrix<double>& A)
{
  // Make sure the matrix is square
  Assert(A.n_block_rows() == A.n_block_cols(),
	 ExcNotImplemented());

  // Define data objects
  std::vector<unsigned int> zero_diag;
  double max_diag(0.0);
  double min_diag(1e+25);
  
  // Find zero element on the diagonal
  for (unsigned int i=0; i < A.m(); i++)
    {
      const double diag_element = A.diag_element(i);
      if (fabs(diag_element) == 0.0)
	zero_diag.push_back(i);
      else if (diag_element > max_diag)
	max_diag = diag_element;
      else if (diag_element < min_diag)
	min_diag = diag_element;
    }
  
  // Reset element in the diagonal
  const double average = (max_diag + min_diag)/2;

  for (unsigned int i=0; i<zero_diag.size(); i++)
    {
      const unsigned int pos = zero_diag[i];
      A.set(pos,pos,average);
    }
}
#ifdef OPENFCST_WITH_PETSC
void
SolverUtils::repair_diagonal(PETScWrappers::MPI::SparseMatrix & A)
{
    // Define data objects
    std::vector<unsigned int> zero_diag;
    double max_diag(0.0);
    double min_diag(1e+25);

    // Find zero element on the diagonal
    std::pair<unsigned int , unsigned int> range = A.local_range();
    for (unsigned int i=range.first; i < range.second; i++)
    {
        const double diag_element = A.diag_element(i);
        if (fabs(diag_element) == 0.0)
            zero_diag.push_back(i);
        else if (diag_element > max_diag)
            max_diag = diag_element;
        else if (diag_element < min_diag)
            min_diag = diag_element;
    }
    // Reset element in the diagonal
    const double average = (max_diag + min_diag)/2;

    for (unsigned int i=0; i<zero_diag.size(); i++)
    {
        const unsigned int pos = zero_diag[i];
        A.set(pos,pos,average);
    }
    A.compress(VectorOperation::insert);
}
#endif



//---------------------------------------------------------------------------
void 
SolverUtils::repair_diagonal(BlockSparseMatrix<double>& A, 
			     FuelCell::ApplicationCore::FEVector& cell_vector,
			     const FuelCell::ApplicationCore::FEVector& solution
			    )
{
  // Make sure the matrix is square
  Assert(A.n_block_rows() == A.n_block_cols(),
	 ExcNotImplemented());

  // Define data objects
  std::vector<unsigned int> zero_diag;
  double max_diag(0.0);
  double min_diag(1e+25);
  
  // Find zero element on the diagonal
  for (unsigned int i=0; i < A.m(); i++)
    {
      const double diag_element = A.diag_element(i);
      if (fabs(diag_element) == 0.0)
	zero_diag.push_back(i);
      else if (diag_element > max_diag)
	max_diag = diag_element;
      else if (diag_element < min_diag)
	min_diag = diag_element;
    }
  
  // Reset element in the diagonal
  const double average = 1; //(max_diag + min_diag)/2;

  for (unsigned int i=0; i<zero_diag.size(); i++)
    {
      const unsigned int pos = zero_diag[i];
      A.set(pos,pos,average);
      cell_vector(pos) = solution(pos); //0;
    }
}
