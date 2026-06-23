//     Copyright (c) 2012 Vadym Kliuchnikov sqct(dot)software(at)gmail(dot)com
//
//     This file is part of SQCT.
//
//     SQCT is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     SQCT is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with SQCT.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef TCOUNT_H
#define TCOUNT_H

#include "rint.h"
#include "matrix2x2.h"

#include <vector>
#include <string>
#include <iostream>

struct min_unitaries
{
  min_unitaries() : min_t_count(-1), k(0), m(0), factor_calls(0), norm_solver_calls(0) {}

  long min_t_count;
  zwt x;
  std::vector<zwt> y;
  int k; // determinant is \w^k
  long m; // power of \sqrt{2} in denominator
  void to_canonical_form();
  bool operator == ( const min_unitaries& rhs ) const;

  static std::string short_title();
  std::string short_str() const;

  long factor_calls;
  long norm_solver_calls;

  operator matrix2x2<mpz_class> () const;


};

std::ostream& operator<< (  std::ostream& out , const min_unitaries& mu );
min_unitaries min_t_count( const zwt& x, long m, int k );


#endif // TCOUNT_H
