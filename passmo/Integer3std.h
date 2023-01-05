/*----------------------------------------------------------------------------*/
 /* Integer3.h                                                                */
 /*                                                                           */
 /* Vecteur d'entiers à 3 dimensions (inspiré de Real3 - Arcane)              */
 /*---------------------------------------------------------------------------*/
 #ifndef INTEGER3_H
 #define INTEGER3_H
 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/

#include <arcane/utils/Real3.h>
#include <arcane/utils/Numeric.h>
#include <arcane/Assertion.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

struct Int3POD {
  public:
	Integer m_i;
	Integer m_j;
	Integer m_k;

   	Int3POD() { m_i = 0; m_j = 0; m_k = 0; }
   	Int3POD(Integer ai,Integer aj,Integer ak) { m_i = ai; m_j = aj; m_k = ak; }
   	Int3POD(const Int3POD& f) { m_i = f.m_i; m_j = f.m_j; m_k = f.m_k; }

   	Int3POD(const Real3& f) { m_i = (Integer)f.x; m_j = (Integer)f.y; m_k = (Integer)f.z; }

   Int3POD& assign(const Int3POD& f)
	{ m_i = f.m_i; m_j = f.m_j; m_k = f.m_k; return (*this); }

	Int3POD&	operator=(Int3POD f){
		return assign(f);
	}

	Integer operator[](Integer i) const {
#ifdef _DEBUG
	assert(i >= 0 && i < 3);// Trying to use an index different than 0, 1 or 2 on a Integer3"
#endif
		return (&m_i)[i];
	}


	Integer& operator[](Integer i) {
#ifdef _DEBUG
	assert(i >= 0 && i < 3);// Trying to use an index different than 0, 1 or 2 on a Integer3"
#endif
		return (&m_i)[i];
	}
 };

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 class Integer3
 : public Int3POD
 {
  public:

   Integer3(): Int3POD() {}
   Integer3(Integer ai,Integer aj,Integer ak):Int3POD(ai,aj,ak) {}
   Integer3(const Integer3& f):Int3POD((const Int3POD&)f) {}
   Integer3(const Real3& f):Int3POD(f) {}

   explicit Integer3(const Int3POD& f):Int3POD(f) {}

   const Integer3&  operator= (Integer v)
   {
     m_i = m_j = m_k = v; return (*this);
   }

   const Integer3& operator=(Real3 f)
   {
     m_i = (Integer)f.x; m_j = (Integer)f.y; m_k = (Integer)f.z;
     return (*this);
   }

  public:

   static Integer3 null() { return {0,0,0}; }
   static Integer3 zero() { return {0,0,0}; }

  public:

   Integer3 copy() const  { return (*this); }
   Integer3& reset() { m_i = m_j = m_k = 0; return (*this); }
   Integer3& assign(Integer ai,Integer aj,Integer ak)
     { m_i = ai; m_j = aj; m_k = ak; return (*this); }
   Integer3& assign(Integer3 f)
     { m_i = f.m_i; m_j = f.m_j; m_k = f.m_k; return (*this); }

   Integer3& operator=(const Integer3& f) { return assign(f); }

   Integer abs2() const
     { return m_i*m_i + m_j*m_j + m_k*m_k; }
   Integer abs() const
     { return sqrt((Real)abs2()); }

   istream& assign(istream& i);
   ostream& print(ostream& o) const;
   ostream& printXyz(ostream& o) const;

   Integer3& add(Integer3 b) { m_i += b.m_i; m_j += b.m_j; m_k += b.m_k; return (*this); }
   Integer3& sub(Integer3 b) { m_i -= b.m_i; m_j -= b.m_j; m_k -= b.m_k; return (*this); }
   Integer3& mul(Integer3 b) { m_i *= b.m_i; m_j *= b.m_j; m_k *= b.m_k; return (*this); }
   Integer3& div(Integer3 b) { m_i/=b.m_i; m_j/=b.m_j; m_k/=b.m_k; return (*this); }
   Integer3& addSame(Integer b) { m_i += b; m_j += b; m_k += b; return (*this); }
   Integer3& subSame(Integer b) { m_i -= b; m_j -= b; m_k -= b; return (*this); }
   Integer3& mulSame(Integer b) { m_i *= b; m_j *= b; m_k *= b; return (*this); }
   Integer3& divSame(Integer b) { m_i /= b; m_j /= b; m_k /= b; return (*this); }
   Integer3& operator+= (Integer3 b) { return add(b); }
   Integer3& operator-= (Integer3 b) { return sub(b); }
   Integer3& operator*= (Integer3 b) { return mul(b); }
   void operator*= (Integer  b) { m_i *= b; m_j *= b; m_k *= b; }
   Integer3& operator/= (Integer3 b) { return div(b); }
   void operator/= (Integer  b) { m_i /= b; m_j /= b; m_k /= b; }
   Integer3 operator+  (Integer3 b)  const { return Integer3(m_i+b.m_i,m_j+b.m_j,m_k+b.m_k); }
   Integer3 operator-  (Integer3 b)  const { return Integer3(m_i-b.m_i,m_j-b.m_j,m_k-b.m_k); }
   Integer3 operator-() const { return Integer3(-m_i,-m_j,-m_k); }
   Integer3 operator*(Integer3 b) const { return Integer3(m_i*b.m_i,m_j*b.m_j,m_k*b.m_k); }
   Integer3 operator/(Integer3 b) const { return Integer3(m_i/b.m_i,m_j/b.m_j,m_k/b.m_k); }

   Integer3& normalize()
   {
     Integer d = abs();
     if (!math::isZero(d))
       divSame(d);
     return (*this);
   }

   bool operator==(Integer3 b) const
   {
     return _eq(m_i,b.m_i) &&  _eq(m_j,b.m_j) && _eq(m_k,b.m_k);
   }
   bool operator!=(Integer3 b) const
     { return !operator==(b); }


  private:

   static bool _eq(Integer a,Integer b)
   {
     return math::isEqual(a,b);
   }
  };

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 inline Integer3 operator*(Integer sca,Integer3 vec) {
   return Integer3(vec.m_i*sca,vec.m_j*sca,vec.m_k*sca);
 }

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 inline Integer3 operator*(Integer3 vec,Integer sca) {
   return Integer3(vec.m_i*sca,vec.m_j*sca,vec.m_k*sca);
 }

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 inline Integer3 operator/(Integer3 vec,Integer sca) {
   return Integer3(vec.m_i/sca,vec.m_j/sca,vec.m_k/sca);
 }

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 inline bool operator<(Integer3 v1,Integer3 v2) {
	 Integer i1 = v1.m_i;
	 Integer j1 = v1.m_j;
	 Integer k1 = v1.m_k;

   if (i1 == v2.m_i) {

	   if (j1 == v2.m_j)

		   return k1 < v2.m_k;
     else

    	 return j1 < v2.m_j;
   }
   return (i1 < v2.m_i);
 }

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/
 inline ostream& operator<< (ostream& o,Integer3 t) {
   return t.printXyz(o);
 }
 inline istream& operator>> (istream& i,Integer3& t) {
   return t.assign(i);
 }

 /*---------------------------------------------------------------------------*/
 /*---------------------------------------------------------------------------*/

 #endif
