//

#include "../include/ghf.h"

using namespace std;
using namespace tensor_hao;

Ghf::Ghf() { }

Ghf::~Ghf() { }

void Ghf::run()
{
    initialParameters();

    for(size_t i = 0; i < ( method.annealStep+1 ); ++i)
    {
        selfConsistentLoop();

        if( variationalEnergy < minimumEnergy )
        {
            minimumEnergy = variationalEnergy;
            minimumState = variationalState;
            cout<<"The minimumEnergy is updated to "<<fixed<<setprecision(16)<<minimumEnergy<<"\n"<<endl;
            writeFile(minimumEnergy, "minimumEnergy.dat");
            minimumState.write("minimumState.dat");
        }
        else
        {
            cout<<"The minimumEnergy is not updated! "<<fixed<<setprecision(16)<<minimumEnergy<<"\n"<<endl;
        }

        if( i!=method.annealStep )
        {
            double anneal = method.annealMagnitude*(method.annealStep-i)/method.annealStep;
            cout<<"Anneal the order parameter: "<<anneal<<"\n"<<endl;
            annealOrderParameter( anneal );
        }
    }

    prepareStop();
}

void Ghf::initialParameters()
{
    if( MPIRank()==0 ) method.read("ghf_param");
    MPIBcast(method);

    randomHaoInit(method.seed, 1);
    if( method.seed != -1 ) randomHaoSave();

    if( MPIRank()==0 ) model.read("model_param");
    MPIBcast(model);

    setH0();

    minimumEnergy = 1e300;

    initialHartreeFockEssential();
}

void Ghf::selfConsistentLoop()
{
    size_t i;
    for(i = 0; i < method.maxIterateStep; ++i)
    {
        oneSelfConsistentStep();
        if( convergeValue < method.convergeTolerance ) break;
    }
    if( i < method.maxIterateStep ) cout<<"\nIt take "<<i<<" step to converge.\n"<<endl;
    else cout<<"\nWARNING!!! Does not converge!!! \n"<<endl;
}

void Ghf::oneSelfConsistentStep()
{
    double variationalEnergyBefore= variationalEnergy;

    setMeanFieldVectorsValuesFromOrderParameter();
    setVariationalStateFromMeanFieldVectors();
    backupAndUpdateOrderParameterFromVariationalState();
    setVariationalEnergyFromOrderParameterAndMeanFieldValues();
    relaxOrderParameter();

    if( method.convergeType == "energy" )
    {
        double difference = abs( variationalEnergy-variationalEnergyBefore );
        convergeValue = difference / abs( variationalEnergy );
    }
    else if( method.convergeType == "orderParameter" )
    {
        size_t L  = model.getL();

        double difference(0), total(0);
        for(size_t i = 0; i < 2*L; ++i)
        {
            for(size_t j = 0; j < 2*L; ++j)
            {
                difference += abs( densityMatrix(j,i) - densityMatrixBefore(j,i) );
                total += abs( densityMatrix(j,i) );
            }

        }
        convergeValue = difference / total;
    }
}

void Ghf::prepareStop()
{
    SDSDOperation sdsdOperation(minimumState, minimumState);
    HubbardKanamoriMeasureObserveSDSD meas(model);
    meas.addMeasurement(sdsdOperation, 1.0);
    meas.write();

    writeNplusOneMinusOneWalker();
}
