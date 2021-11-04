//

#include "../include/ghf.h"

using namespace std;
using namespace tensor_hao;

void Ghf::setH0()
{
    H0 = model.getT();
}

void Ghf::initialHartreeFockEssential()
{
    size_t L = model.getL();
    densityMatrixBefore.resize( 2*L, 2*L );

    if( method.initialType == "setFromModel" )
    {
        densityMatrix.resize(2*L, 2*L);
        densityMatrix = complex<double>(0,0);
    }
    else if( method.initialType == "readWaveFunction" )
    {
        variationalState.read("phi.dat");

        HubbardKanamoriMeasureFixSDSD meas(model, variationalState);

        SDSDOperation sdsdOperation(variationalState, variationalState);

        meas.addMeasurement(sdsdOperation, 1.0);

        minimumEnergy = ( meas.returnEnergy() ).real();
        minimumState = variationalState;
        cout<<"Read wave function, variational energy: "<<fixed<<setprecision(16)<<minimumEnergy<<endl;

        densityMatrix = sdsdOperation.returnGreenMatrix();
    }
    else if( method.initialType == "readOrderParameter" )
    {
        densityMatrix.read("densityMatrix.dat");
    }
    else
    {
        cout<<"Error!!! Do not know initialType "<<method.initialType<<endl;
        exit(1);
    }

    setMeanFieldVectorsValuesFromOrderParameter();
    setVariationalStateFromMeanFieldVectors();
    backupAndUpdateOrderParameterFromVariationalState();
    setVariationalEnergyFromOrderParameterAndMeanFieldValues();
}

void Ghf::setMeanFieldVectorsValuesFromOrderParameter()
{
    size_t L = model.getL();
    size_t numberOfKana = model.getNumberOfKana();
    const TensorHao<int, 1> &site_i = model.getSite_i();
    const TensorHao<int, 1> &site_j = model.getSite_j();
    const TensorHao< double, 1> &U  = model.getU();
    const TensorHao< double, 1> &U1 = model.getU1();
    const TensorHao< double, 1> &U2 = model.getU2();
    const TensorHao< double, 1> &J  = model.getJ();

    int ki, kj;

    meanFieldVectors = H0;

    //U
    for(size_t i = 0; i < L; ++i)
    {
        meanFieldVectors(i,  i  ) += U(i) * densityMatrix(i+L, i+L);
        meanFieldVectors(i+L,i+L) += U(i) * densityMatrix(i,   i  );
        meanFieldVectors(i,  i+L) -= U(i) * densityMatrix(i+L, i  );
        meanFieldVectors(i+L,i  ) -= U(i) * densityMatrix(i,   i+L);
    }

    //U1
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        meanFieldVectors(ki,   ki  ) += U1(i) * densityMatrix(kj+L, kj+L);
        meanFieldVectors(kj+L, kj+L) += U1(i) * densityMatrix(ki,   ki  );
        meanFieldVectors(kj+L, ki  ) -= U1(i) * densityMatrix(ki,   kj+L);
        meanFieldVectors(ki,   kj+L) -= U1(i) * densityMatrix(kj+L, ki  );

        meanFieldVectors(kj,   kj  ) += U1(i) * densityMatrix(ki+L, ki+L);
        meanFieldVectors(ki+L, ki+L) += U1(i) * densityMatrix(kj,   kj  );
        meanFieldVectors(ki+L, kj  ) -= U1(i) * densityMatrix(kj,   ki+L);
        meanFieldVectors(kj,   ki+L) -= U1(i) * densityMatrix(ki+L, kj  );
    }

    //U2
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        meanFieldVectors(ki, ki) += U2(i) * densityMatrix(kj, kj);
        meanFieldVectors(kj, kj) += U2(i) * densityMatrix(ki, ki);
        meanFieldVectors(ki, kj) -= U2(i) * densityMatrix(kj, ki);
        meanFieldVectors(kj, ki) -= U2(i) * densityMatrix(ki, kj);

        meanFieldVectors(ki+L, ki+L) += U2(i) * densityMatrix(kj+L, kj+L);
        meanFieldVectors(kj+L, kj+L) += U2(i) * densityMatrix(ki+L, ki+L);
        meanFieldVectors(ki+L, kj+L) -= U2(i) * densityMatrix(kj+L, ki+L);
        meanFieldVectors(kj+L, ki+L) -= U2(i) * densityMatrix(ki+L, kj+L);
    }

    //J
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        meanFieldVectors(ki,   kj  ) += J(i) * densityMatrix(kj+L, ki+L);
        meanFieldVectors(kj+L, ki+L) += J(i) * densityMatrix(ki,   kj  );
        meanFieldVectors(ki,   ki+L) -= J(i) * densityMatrix(kj+L, kj  );
        meanFieldVectors(kj+L, kj  ) -= J(i) * densityMatrix(ki,   ki+L);

        meanFieldVectors(ki,   kj  ) += J(i) * densityMatrix(ki+L, kj+L);
        meanFieldVectors(ki+L, kj+L) += J(i) * densityMatrix(ki,   kj  );
        meanFieldVectors(ki,   kj+L) -= J(i) * densityMatrix(ki+L, kj  );
        meanFieldVectors(ki+L, kj  ) -= J(i) * densityMatrix(ki,   kj+L);

        meanFieldVectors(kj,   ki  ) += J(i) * densityMatrix(ki+L, kj+L);
        meanFieldVectors(ki+L, kj+L) += J(i) * densityMatrix(kj,   ki  );
        meanFieldVectors(kj,   kj+L) -= J(i) * densityMatrix(ki+L, ki  );
        meanFieldVectors(ki+L, ki  ) -= J(i) * densityMatrix(kj,   kj+L);

        meanFieldVectors(kj,   ki  ) += J(i) * densityMatrix(kj+L, ki+L);
        meanFieldVectors(kj+L, ki+L) += J(i) * densityMatrix(kj,   ki  );
        meanFieldVectors(kj,   ki+L) -= J(i) * densityMatrix(kj+L, ki  );
        meanFieldVectors(kj+L, ki  ) -= J(i) * densityMatrix(kj,   ki+L);

    }

    meanFieldValues.resize(2*L);

    BL_NAME(eigen)(meanFieldVectors, meanFieldValues);
}

void Ghf::setVariationalStateFromMeanFieldVectors()
{
    size_t L = model.getL();
    size_t N = model.getN();

    variationalState.resize(2*L, N);
    TensorHao<complex<double>,2> &wf = variationalState.wfRef();
    copy( meanFieldVectors.data(), meanFieldVectors.data()+2*L*N, wf.data() );
}

void Ghf::backupAndUpdateOrderParameterFromVariationalState()
{
    densityMatrixBefore = move(densityMatrix);

    SDSDOperation sdsdOperation(variationalState, variationalState);
    densityMatrix = sdsdOperation.returnGreenMatrix();
}

void Ghf::setVariationalEnergyFromOrderParameterAndMeanFieldValues()
{
    variationalEnergy = 0.0;

    size_t L = model.getL();
    size_t N = model.getN();
    size_t numberOfKana = model.getNumberOfKana();
    const TensorHao<int, 1> &site_i = model.getSite_i();
    const TensorHao<int, 1> &site_j = model.getSite_j();
    const TensorHao< double, 1> &U  = model.getU();
    const TensorHao< double, 1> &U1 = model.getU1();
    const TensorHao< double, 1> &U2 = model.getU2();
    const TensorHao< double, 1> &J  = model.getJ();

    int ki, kj;

    for(size_t i = 0; i < N; ++i) variationalEnergy += meanFieldValues(i);

    for(size_t i = 0; i < L; ++i)
    {
        variationalEnergy += U(i) * ( -densityMatrixBefore(i,   i  )*densityMatrix(i+L, i+L)
                                      -densityMatrixBefore(i+L, i+L)*densityMatrix(i,   i  )
                                      +densityMatrixBefore(i,   i+L)*densityMatrix(i+L, i  )
                                      +densityMatrixBefore(i+L, i  )*densityMatrix(i,   i+L)
                                      +densityMatrix(i, i  )*densityMatrix(i+L, i+L)
                                      -densityMatrix(i, i+L)*densityMatrix(i+L, i )
                                    ).real();
    }

    //U1
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        variationalEnergy += U1(i) * ( -densityMatrixBefore(ki,   ki  ) * densityMatrix(kj+L, kj+L)
                                       -densityMatrixBefore(kj+L, kj+L) * densityMatrix(ki,   ki  )
                                       +densityMatrixBefore(kj+L, ki  ) * densityMatrix(ki,   kj+L)
                                       +densityMatrixBefore(ki,   kj+L) * densityMatrix(kj+L, ki  )
                                       -densityMatrixBefore(kj,   kj  ) * densityMatrix(ki+L, ki+L)
                                       -densityMatrixBefore(ki+L, ki+L) * densityMatrix(kj,   kj  )
                                       +densityMatrixBefore(ki+L, kj  ) * densityMatrix(kj,   ki+L)
                                       +densityMatrixBefore(kj,   ki+L) * densityMatrix(ki+L, kj  )
                                       +densityMatrix(ki,  ki  )*densityMatrix(kj+L,kj+L)
                                       -densityMatrix(ki,  kj+L)*densityMatrix(kj+L,ki  )
                                       +densityMatrix(ki+L,ki+L)*densityMatrix(kj  ,kj  )
                                       -densityMatrix(ki+L,kj  )*densityMatrix(kj  ,ki+L)
                                     ).real();
    }

    //U2
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        variationalEnergy += U2(i) * (-densityMatrixBefore(ki,   ki  ) * densityMatrix(kj,   kj  )
                                      -densityMatrixBefore(kj,   kj  ) * densityMatrix(ki,   ki  )
                                      +densityMatrixBefore(ki,   kj  ) * densityMatrix(kj,   ki  )
                                      +densityMatrixBefore(kj,   ki  ) * densityMatrix(ki,   kj  )
                                      -densityMatrixBefore(ki+L, ki+L) * densityMatrix(kj+L, kj+L)
                                      -densityMatrixBefore(kj+L, kj+L) * densityMatrix(ki+L, ki+L)
                                      +densityMatrixBefore(ki+L, kj+L) * densityMatrix(kj+L, ki+L)
                                      +densityMatrixBefore(kj+L, ki+L) * densityMatrix(ki+L, kj+L)
                                      +densityMatrix(ki,  ki  )*densityMatrix(kj,  kj  )
                                      -densityMatrix(ki,  kj  )*densityMatrix(kj,  ki  )
                                      +densityMatrix(ki+L,ki+L)*densityMatrix(kj+L,kj+L)
                                      -densityMatrix(ki+L,kj+L)*densityMatrix(kj+L,ki+L)
                                     ).real();
    }

    //J
    for(size_t i = 0; i < numberOfKana; ++i)
    {
        ki = site_i(i); kj = site_j(i);

        variationalEnergy += J(i) * (- densityMatrixBefore(ki,   kj  ) * densityMatrix(kj+L, ki+L)
                                     - densityMatrixBefore(kj+L, ki+L) * densityMatrix(ki,   kj  )
                                     + densityMatrixBefore(ki,   ki+L) * densityMatrix(kj+L, kj  )
                                     + densityMatrixBefore(kj+L, kj  ) * densityMatrix(ki,   ki+L)
                                     - densityMatrixBefore(ki,   kj  ) * densityMatrix(ki+L, kj+L)
                                     - densityMatrixBefore(ki+L, kj+L) * densityMatrix(ki,   kj  )
                                     + densityMatrixBefore(ki,   kj+L) * densityMatrix(ki+L, kj  )
                                     + densityMatrixBefore(ki+L, kj  ) * densityMatrix(ki,   kj+L)
                                     - densityMatrixBefore(kj,   ki  ) * densityMatrix(ki+L, kj+L)
                                     - densityMatrixBefore(ki+L, kj+L) * densityMatrix(kj,   ki  )
                                     + densityMatrixBefore(kj,   kj+L) * densityMatrix(ki+L, ki  )
                                     + densityMatrixBefore(ki+L, ki  ) * densityMatrix(kj,   kj+L)
                                     - densityMatrixBefore(kj,   ki  ) * densityMatrix(kj+L, ki+L)
                                     - densityMatrixBefore(kj+L, ki+L) * densityMatrix(kj,   ki  )
                                     + densityMatrixBefore(kj,   ki+L) * densityMatrix(kj+L, ki  )
                                     + densityMatrixBefore(kj+L, ki  ) * densityMatrix(kj,   ki+L)
                                     +densityMatrix(ki,kj  )*densityMatrix(kj+L,ki+L)
                                     -densityMatrix(ki,ki+L)*densityMatrix(kj+L,kj  )
                                     +densityMatrix(ki,kj  )*densityMatrix(ki+L,kj+L)
                                     -densityMatrix(ki,kj+L)*densityMatrix(ki+L,kj  )
                                     +densityMatrix(kj,ki  )*densityMatrix(ki+L,kj+L)
                                     -densityMatrix(kj,kj+L)*densityMatrix(ki+L,ki  )
                                     +densityMatrix(kj,ki  )*densityMatrix(kj+L,ki+L)
                                     -densityMatrix(kj,ki+L)*densityMatrix(kj+L,ki  )
                                    ).real();
    }

    cout<<"Variational energy: "<<fixed<<setprecision(16)<<variationalEnergy<<endl;

//    Test code
//    HubbardKanamoriMeasureFixSDSD meas(model, variationalState);
//    SDSDOperation sdsdOperation(variationalState, variationalState);
//    meas.addMeasurement(sdsdOperation, 1.0);
//    double mE = ( meas.returnEnergy() ).real();
//    if( abs( mE-variationalEnergy )> 1e-12 )
//    {
//        cout<<"Error!!!"<<endl;
//        cout<<"Variational energy: "<<fixed<<setprecision(16)<<mE<<endl;
//    }
//    cout<<endl;
}

void Ghf::relaxOrderParameter()
{
    size_t L = model.getL();
    for(size_t i = 0; i < 2*L; ++i)
    {
        for(size_t j = 0; j < 2*L; ++j)
        {
            densityMatrix(j,i) = densityMatrixBefore(j,i)
                                 + method.relaxMagnitude * ( densityMatrix(j,i)-densityMatrixBefore(j,i) );

        }
    }
}

void Ghf::annealOrderParameter(double anneal)
{
    size_t L = model.getL();
    for(size_t i = 0; i < 2*L; ++i)
    {
        for(size_t j = i; j < 2*L; ++j)
        {
            if( abs( densityMatrix(j,i)  )> 1e-12 )
            {
                densityMatrix(j,i) *= ( 1.0 + (2.0*uniformHao()-1.0 )*anneal );
                densityMatrix(i,j) = conj( densityMatrix(j,i) );
            }
            else
            {
                densityMatrix(j,i) = 0.0;
                densityMatrix(i,j) = 0.0;
            }
        }
    }
}

void Ghf::writeNplusOneMinusOneWalker()
{
    size_t L = model.getL();
    size_t N = model.getN();

    SD plusOneWf(2*L, N+1), minusOneWf(2*L, N-1) ;

    TensorHao<complex<double>,2> &wfplus = plusOneWf.wfRef();
    copy( meanFieldVectors.data(), meanFieldVectors.data()+2*L*(N+1), wfplus.data() );
    plusOneWf.write("phiPlus.dat");

    TensorHao<complex<double>,2> &wfminus = minusOneWf.wfRef();
    copy( meanFieldVectors.data(), meanFieldVectors.data()+2*L*(N-1), wfminus.data() );
    minusOneWf.write("phiMinus.dat");
}
