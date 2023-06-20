#pragma once

#include "xEval.h"
#include "xVector.h"

#include "xGenericSparseMatrix.h"


struct xGetRealPart : public std::unary_function< std::complex<double> , double >
{
public:
    typedef std::unary_function< std::complex<double>, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return t.real(); }
    //    const result_type& operator()( const argument_type& t ) const { return const_cast<double&>(t.real()); }
};


struct xGetImagPart : public std::unary_function< std::complex<double> , double >
{
public:
    typedef std::unary_function< std::complex<double>, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return t.imag(); }
};



struct xGetComplex : public std::unary_function< std::complex<double> , std::complex<double> >
{
public:
    typedef std::unary_function< std::complex<double>, std::complex<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return t; }
    //const result_type& operator()( const argument_type& t ) const { return t; }
};

struct xGetPressureIndB : public std::unary_function< std::complex<double> , double >
{
public:
    typedef std::unary_function< std::complex<double>, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return 20.*log10(abs(t)/0.00002); }
    //const result_type& operator()( const argument_type& t ) const { return t; }
};


struct xReflection : public std::unary_function< std::complex<double> , std::complex<double> >
{
public:
    typedef std::unary_function< std::complex<double>, std::complex<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return (t-415.1)/(t+415.1); }
    //const result_type& operator()( const argument_type& t ) const { return t; }
};



struct xGetNormalvelocity : public std::unary_function< xtensor::xVectorDoubleComplex, std::complex<double> >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, std::complex<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;
    

    xGetNormalvelocity(double omega_) : omega(omega_) {};
    result_type operator()( argument_type& t) const { return (-t[0]/(1.213*j*omega)); }
    private:
    const std::complex<double> j{0.,1.};
    double omega;
};


template <typename T>
    struct xGetModu : public std::unary_function< T, double >
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return std::abs(t); }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    };


template <typename T>
    struct xGetModuV : public std::unary_function< T, double >
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return sqrt(pow(std::abs(t[0]/(j*18400.)),2)); }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    private:
    const std::complex<double> j{0.,1.};
    
    };

template <typename T>
    struct xGetPhase : public std::unary_function< T, double >
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return std::arg(t); }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    };


template <typename T>
    struct xGetSquare : public std::unary_function<T, double>
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return t*t; }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    };


 template <typename T>
    struct xAbsorption : public std::unary_function<T, double>
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return 1-std::abs(t)*std::abs(t); }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    };


 template <typename T>
    struct xGetSqrt : public std::unary_function<T, double>
    {
    public:
        typedef std::unary_function< T, double >               base_type;
        typedef typename base_type::argument_type         argument_type;
        typedef typename base_type::result_type           result_type;

        result_type operator()( argument_type& t ) const { return sqrt(t); }
        //const result_type& operator()( const argument_type& t ) const { return t; }
    };


template<class T1, class T2, class T3>
    class xMinusComplex : public std::binary_function< T1, T2, T3 >
    {
    public:
        typedef typename std::binary_function<T1, T2, T3 >::result_type result_type;
        typedef typename std::binary_function<T1, T2, T3 >::first_argument_type first_argument_type;
        typedef typename std::binary_function<T1, T2, T3 >::second_argument_type second_argument_type;

        result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f-s; }
    };

template<class T1, class T2, class T3>
    class xMinusComplexV : public std::binary_function< T1, T2, T3 >
    {
    public:
        typedef typename std::binary_function<T1, T2, T3 >::result_type result_type;
        typedef typename std::binary_function<T1, T2, T3 >::first_argument_type first_argument_type;
        typedef typename std::binary_function<T1, T2, T3 >::second_argument_type second_argument_type;

        result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return xtensor::xVector<>(f[0]-s[0], f[1]-s[1]); }
    };

template<class T1, class T2, class T3>
    class xPlusComplexVS : public std::binary_function< T1, T2, T3 >
    {
    public:
        typedef typename std::binary_function<T1, T2, T3 >::result_type result_type;
        typedef typename std::binary_function<T1, T2, T3 >::first_argument_type first_argument_type;
        typedef typename std::binary_function<T1, T2, T3 >::second_argument_type second_argument_type;

        result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f+s[0]+s[1]; }
    };


template<class T1, class T2, class T3>
    class xDivisComplex : public std::binary_function< T1, T2, T3 >
    {
    public:
        typedef typename std::binary_function<T1, T2, T3 >::result_type result_type;
        typedef typename std::binary_function<T1, T2, T3 >::first_argument_type first_argument_type;
        typedef typename std::binary_function<T1, T2, T3 >::second_argument_type second_argument_type;

        result_type operator()(const first_argument_type& f, const second_argument_type& s) const {return f/s; }
    };



struct xGetRealPartV : public std::unary_function< xtensor::xVectorDoubleComplex , xtensor::xVector<double> >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, xtensor::xVector<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return xtensor::xVector<>(t[0].real(), t[1].real(), t[2].real()); }
};

struct xGetConjugV : public std::unary_function< xtensor::xVectorDoubleComplex , xtensor::xVectorDoubleComplex >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, xtensor::xVectorDoubleComplex >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return xtensor::xVector<std::complex<double>>(std::conj(t[0]), std::conj(t[1]), std::conj(t[2])); }
};


struct xGetNormTensor : public std::unary_function< xtensor::xTensor2<std::complex<double> > , double >
{
public:
    typedef std::unary_function< xtensor::xTensor2<std::complex<double> >, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return t(0,0).real(); }
};



struct xGetRealPartVComponent : public std::unary_function< xtensor::xVectorDoubleComplex , double >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    xGetRealPartVComponent(int comp) : component(comp) {};

    result_type operator()( argument_type& t ) const { return t[component].real(); }

    private:
    const int component;
};

struct xGetNormVComponent : public std::unary_function< xtensor::xVectorDoubleComplex , double >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, double >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    xGetNormVComponent(int comp) : component(comp) {};

    result_type operator()( argument_type& t ) const { return std::abs(t[component]); }

    private:
    const int component;
};


struct xGetImagPartV : public std::unary_function< xtensor::xVectorDoubleComplex , xtensor::xVector<double> >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, xtensor::xVector<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return xtensor::xVector<>(t[0].imag(), t[1].imag(), t[2].imag()); }
};

struct xGetComplexV : public std::unary_function< xtensor::xVectorDoubleComplex , std::complex<double> >
{
public:
    typedef std::unary_function< xtensor::xVectorDoubleComplex, std::complex<double> >               base_type;
    typedef typename base_type::argument_type         argument_type;
    typedef typename base_type::result_type           result_type;

    result_type operator()( argument_type& t ) const { return t[0]; }
};

template<class UnaryOperator>
class xEvalVelocity : public xfem::xEval<typename UnaryOperator::result_type> {
public :
    typedef typename UnaryOperator::result_type   result_type;
    xEvalVelocity(xfem::xEval<xtensor::xVector<std::complex<double> > > &eval_grad_,
                  double rho_, double omega_) : eval_grad(eval_grad_), rho(rho_), omega(omega_), funct() {}

    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, result_type& res) const {

        xtensor::xVector<std::complex<double> > grad;
        eval_grad(geo_appro, geo_integ, grad);
        xtensor::xVector<std::complex<double> > velocity = j / (rho * omega) * grad;

        res = funct(velocity);
    }
    

private:
    xfem::xEval<xtensor::xVector<std::complex<double> > > &eval_grad;
    double rho, omega;
    const std::complex<double> j{0.,1.};
    const UnaryOperator funct;
};


template<class UnaryOperator>
class xEvalDisp : public xfem::xEval<typename UnaryOperator::result_type> {
public :
    typedef typename UnaryOperator::result_type   result_type;
    xEvalDisp(xfem::xEval<xtensor::xVector<std::complex<double> > > &eval_field_) : eval_field(eval_field_), funct() {}

    void operator()(const xfem::xGeomElem*  geo_appro, const xfem::xGeomElem* geo_integ, result_type& res) const {

        xtensor::xVector<std::complex<double> > field;
        eval_field(geo_appro, geo_integ, field);
        xtensor::xVector<std::complex<double> > disp = field;

        res = funct(disp);
    }
    

private:
    xfem::xEval<xtensor::xVector<std::complex<double> > > &eval_field;
    const std::complex<double> j{0.,1.};
    const UnaryOperator funct;
};



// TO BE REMOVED WHEN THE COMPLEX VERSION OF MUMPS WILL BE AVAILABLE ON TITAN !!!
template < class MT, class TRIP >
void fillTripletSym(MT & AA, std::vector<TRIP > &Triplets)
{
    int k,l,ke,le;
    ke = xlinalg::xPolicyGenericSparseMatrix < typename MT::matrix_pattern, typename MT::matrix_storage, typename MT::matrix_indexing >::firstIndexEnd(AA.getN(),AA.getM(), AA.getNNZ(),AA.index1,AA.index2);
    std::cout<<"KE IS "<<ke<<std::endl;
    for (k = 0; k < ke; ++k)
    {
        l = xlinalg::xPolicyGenericSparseMatrix < typename MT::matrix_pattern, typename MT::matrix_storage, typename MT::matrix_indexing  >::secondIndexStart(k);
        le = xlinalg::xPolicyGenericSparseMatrix < typename MT::matrix_pattern, typename MT::matrix_storage, typename MT::matrix_indexing  >::secondIndexEnd(k);
        for (; l < le; ++l)
        {

            auto ii = xlinalg::xPolicyGenericSparseMatrix < typename MT::matrix_pattern, typename MT::matrix_storage, typename MT::matrix_indexing >::indexI(k,l)-1;
            auto jj = xlinalg::xPolicyGenericSparseMatrix < typename MT::matrix_pattern, typename MT::matrix_storage, typename MT::matrix_indexing >::indexJ(k,l)-1;

            //            cout<<ii<<" "<<jj<<endl;

            Triplets.push_back(TRIP(ii,jj, AA(ii+1, jj+1)));

            if(ii != jj) Triplets.push_back(TRIP(jj,ii, AA(ii+1, jj+1)));

            //            os<<xPolicyGenericSparseMatrix < PATTERN,STORAGE_TYPE,INDEXING >::indexI(k,l);
            //            os<<" "<<xPolicyGenericSparseMatrix < PATTERN,STORAGE_TYPE,INDEXING >::indexJ(k,l);
            //            os<<" "<<( (T * ) data )[xPolicyGenericSparseMatrix < PATTERN,STORAGE_TYPE,INDEXING > ::indexVal(k,l)]<<std::endl;
        }
    }

}
