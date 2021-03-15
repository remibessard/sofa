/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once
#include <SofaBaseLinearSolver/config.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/BaseMatrix.h>

namespace sofa::component::linearsolver
{

	template<int TN> class bloc_index_func
	{
	public:
		enum { N = TN };
		static void split(int& index, int& modulo)
		{
			modulo = index % N;
			index = index / N;
		}
	};

	template<> class bloc_index_func<1>
	{
	public:
		enum { N = 1 };
		static void split(int&, int&)
		{
		}
	};

	template<> class bloc_index_func<2>
	{
	public:
		enum { N = 2 };
		static void split(int& index, int& modulo)
		{
			modulo = index & 1;
			index = index >> 1;
		}
	};

	template<> class bloc_index_func<4>
	{
	public:
		enum { N = 2 };
		static void split(int& index, int& modulo)
		{
			modulo = index & 3;
			index = index >> 2;
		}
	};

	template<> class bloc_index_func<8>
	{
	public:
		enum { N = 2 };
		static void split(int& index, int& modulo)
		{
			modulo = index & 7;
			index = index >> 3;
		}
	};


	// by default, supposing T is a defaulttype::Mat (usefull for type derivated from defaulttype::Mat)
	template<class T>
	class matrix_bloc_traits
	{
	public:
		typedef T Bloc;
		typedef T BlocTranspose;

		typedef typename T::Real Real;
		enum { NL = T::nbLines };
		enum { NC = T::nbCols };
		static const Real& v(const Bloc& b, int row, int col) { return b[row][col]; }
		static void vset(Bloc& b, int row, int col, Real val) { b[row][col] = val; }
		static void vadd(Bloc& b, int row, int col, Real val) { b[row][col] += val; }
		static void clear(Bloc& b) { b.clear(); }
		static bool empty(const Bloc& b)
		{
			for (int i = 0; i < NL; ++i)
				for (int j = 0; j < NC; ++j)
					if (b[i][j] != 0) return false;
			return true;
		}
		static void invert(Bloc& result, const Bloc& b) { result.invert(b); }

		static BlocTranspose transposed(const Bloc& b) { return b.transposed(); }

		static void transpose(BlocTranspose& res, const Bloc& b) { res.transpose(b); }

		static void split_row_index(int& index, int& modulo) { bloc_index_func<NL>::split(index, modulo); }
		static void split_col_index(int& index, int& modulo) { bloc_index_func<NC>::split(index, modulo); }

		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return matrix_bloc_traits<Real>::getElementType(); }
		static const char* Name();
	};

	template <int L, int C, class real>
	class matrix_bloc_traits < defaulttype::Mat<L, C, real> >
	{
	public:
		typedef defaulttype::Mat<L, C, real> Bloc;
		typedef defaulttype::Mat<C, L, real> BlocTranspose;

		typedef real Real;
		enum { NL = L };
		enum { NC = C };
		static const Real& v(const Bloc& b, int row, int col) { return b[row][col]; }
		static void vset(Bloc& b, int row, int col, Real val) { b[row][col] = val; }
		static void vadd(Bloc& b, int row, int col, Real val) { b[row][col] += val; }
		static void clear(Bloc& b) { b.clear(); }
		static bool empty(const Bloc& b)
		{
			for (int i = 0; i < NL; ++i)
				for (int j = 0; j < NC; ++j)
					if (b[i][j] != 0) return false;
			return true;
		}
		static void invert(Bloc& result, const Bloc& b) { result.invert(b); }

		static BlocTranspose transposed(const Bloc& b) { return b.transposed(); }

		static void transpose(BlocTranspose& res, const Bloc& b) { res.transpose(b); }

		static void split_row_index(int& index, int& modulo) { bloc_index_func<NL>::split(index, modulo); }
		static void split_col_index(int& index, int& modulo) { bloc_index_func<NC>::split(index, modulo); }

		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return matrix_bloc_traits<Real>::getElementType(); }
		static const char* Name();
	};

	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<1, 1, float > >::Name() { return "1f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<1, 1, double> >::Name() { return "1d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<2, 2, float > >::Name() { return "2f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<2, 2, double> >::Name() { return "2d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<2, 3, float > >::Name() { return "2x3f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<2, 3, double> >::Name() { return "2x3d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 2, float > >::Name() { return "3x2f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 2, double> >::Name() { return "3x2d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 3, float > >::Name() { return "3f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 3, double> >::Name() { return "3d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<4, 4, float > >::Name() { return "4f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<4, 4, double> >::Name() { return "4d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<6, 6, float > >::Name() { return "6f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<6, 6, double> >::Name() { return "6d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 6, float > >::Name() { return "3x6f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<3, 6, double> >::Name() { return "3x6d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<6, 3, float > >::Name() { return "6x3f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<6, 3, double> >::Name() { return "6x3d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<8, 8, float > >::Name() { return "8f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<8, 8, double> >::Name() { return "8d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<9, 9, float > >::Name() { return "9f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<9, 9, double> >::Name() { return "9d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<12, 12, float > >::Name() { return "12f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Mat<12, 12, double> >::Name() { return "12d"; }

	template <>
	class matrix_bloc_traits < float >
	{
	public:
		typedef float Bloc;
		typedef float Real;
		typedef Bloc BlocTranspose;

		enum { NL = 1 };
		enum { NC = 1 };
		static const Real& v(const Bloc& b, int, int) { return b; }
		static void vset(Bloc& b, int, int, Real val) { b = val; }
		static void vadd(Bloc& b, int, int, Real val) { b += val; }
		static void clear(Bloc& b) { b = 0; }
		static bool empty(const Bloc& b)
		{
			return b == 0;
		}
		static void invert(Bloc& result, const Bloc& b) { result = 1.0f / b; }

		static Bloc transposed(const Bloc& b) { return b; }

		static void transpose(Bloc& res, const Bloc& b) { res = b; }

		static void split_row_index(int& index, int& modulo) { bloc_index_func<NL>::split(index, modulo); }
		static void split_col_index(int& index, int& modulo) { bloc_index_func<NC>::split(index, modulo); }

		static const char* Name() { return "f"; }
		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return sofa::defaulttype::BaseMatrix::ELEMENT_FLOAT; }
		static std::size_t getElementSize() { return sizeof(Real); }
	};

	template <>
	class matrix_bloc_traits < double >
	{
	public:
		typedef double Bloc;
		typedef double Real;
		typedef Bloc BlocTranspose;

		enum { NL = 1 };
		enum { NC = 1 };
		static const Real& v(const Bloc& b, int, int) { return b; }
		static void vset(Bloc& b, int, int, Real val) { b = val; }
		static void vadd(Bloc& b, int, int, Real val) { b += val; }
		static void clear(Bloc& b) { b = 0; }
		static bool empty(const Bloc& b)
		{
			return b == 0;
		}
		static void invert(Bloc& result, const Bloc& b) { result = 1.0 / b; }

		static Bloc transposed(const Bloc& b) { return b; }

		static void transpose(Bloc& res, const Bloc& b) { res = b; }

		static void split_row_index(int& index, int& modulo) { bloc_index_func<NL>::split(index, modulo); }
		static void split_col_index(int& index, int& modulo) { bloc_index_func<NC>::split(index, modulo); }

		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return sofa::defaulttype::BaseMatrix::ELEMENT_FLOAT; }
		static const char* Name() { return "d"; }
	};

	template <>
	class matrix_bloc_traits < int >
	{
	public:
		typedef float Bloc;
		typedef float Real;
		typedef Bloc BlocTranspose;

		enum { NL = 1 };
		enum { NC = 1 };
		static const Real& v(const Bloc& b, int, int) { return b; }
		static void vset(Bloc& b, int, int, Real val) { b = val; }
		static void vadd(Bloc& b, int, int, Real val) { b += val; }
		static void clear(Bloc& b) { b = 0; }
		static bool empty(const Bloc& b)
		{
			return b == 0;
		}
		static void invert(Bloc& result, const Bloc& b) { result = 1.0f / b; }

		static Bloc transposed(const Bloc& b) { return b; }

		static void transpose(Bloc& res, const Bloc& b) { res = b; }

		static void split_row_index(int& index, int& modulo) { bloc_index_func<NL>::split(index, modulo); }
		static void split_col_index(int& index, int& modulo) { bloc_index_func<NC>::split(index, modulo); }

		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return sofa::defaulttype::BaseMatrix::ELEMENT_INT; }
		static const char* Name() { return "f"; }
	};

	template <int N, class T>
	class matrix_bloc_traits < sofa::defaulttype::Vec<N, T> >
	{
	public:
		typedef sofa::defaulttype::Vec<N, T> Bloc;
		typedef T Real;
		typedef Bloc BlocTranspose;

		enum { NL = 1 };
		enum { NC = N };

		static const Real& v(const Bloc& b, int /*row*/, int col) { return b[col]; }
		static void vset(Bloc& b, int /*row*/, int col, Real v) { b[col] = v; }
		static void vadd(Bloc& b, int /*row*/, int col, Real v) { b[col] += v; }
		static void clear(Bloc& b) { b.clear(); }
		static bool empty(const Bloc& b)
		{
			for (int i = 0; i < NC; ++i)
				if (b[i] != 0) return false;
			return true;
		}

		static Bloc transposed(const Bloc& b) { return b; }

		static void transpose(Bloc& res, const Bloc& b) { res = b; }

		static sofa::defaulttype::BaseMatrix::ElementType getElementType() { return matrix_bloc_traits<Real>::getElementType(); }
		static const char* Name();
	};

	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<1, float > >::Name() { return "V1f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<1, double> >::Name() { return "V1d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<2, float > >::Name() { return "V2f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<2, double> >::Name() { return "V2d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<3, float > >::Name() { return "V3f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<3, double> >::Name() { return "V3d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<4, float > >::Name() { return "V4f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<4, double> >::Name() { return "V4d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<5, float > >::Name() { return "V5f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<5, double> >::Name() { return "V5d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<6, float > >::Name() { return "V6f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<6, double> >::Name() { return "V6d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<7, float > >::Name() { return "V7f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<7, double> >::Name() { return "V7d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<8, float > >::Name() { return "V8f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<8, double> >::Name() { return "V8d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<9, float > >::Name() { return "V9f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<9, double> >::Name() { return "V9d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<10, float > >::Name() { return "V10f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<10, double> >::Name() { return "V10d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<11, float > >::Name() { return "V11f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<11, double> >::Name() { return "V11d"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<12, float > >::Name() { return "V12f"; }
	template<> inline const char* matrix_bloc_traits<defaulttype::Vec<12, double> >::Name() { return "V12d"; }

} // namespace sofa::component::linearsolver
