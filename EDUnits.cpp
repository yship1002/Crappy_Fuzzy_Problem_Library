#include "EDUnits.h"
#include <vector>
#include <cmath>


// NEEDED: no division, simplier formation, axu var
EDUnits::EDUnits(BranchingStrategy branchingStrategy) : STModel() {
    this->branching_strategy = branchingStrategy;
    this->scenario_names = {
        ScenarioNames::SCENARIO1,
        ScenarioNames::SCENARIO2,
        ScenarioNames::SCENARIO3
    };
    this->scenario_name = ScenarioNames::SCENARIO1;
    this->probability = 1.0;

    // this->first_stage_IX = {
    //     mc::Interval(0.2229377820190313, 0.2229377820190313),      // molNStream1 0
    //     mc::Interval(123.68485034381521,123.68485034381521),     // molWStream1 1
    //     mc::Interval(0.06147995475275121,0.06147995475275121),      // molNStream2 2
    //     mc::Interval(119.50012464909788,119.50012464909788),     // molWStream2 3
    //     mc::Interval(0.16145782726628008,0.16145782726628008),      // molNStream3 4
    //     mc::Interval(4.184725694717343,4.184725694717343),     // molWStream3 5
    //     mc::Interval(9.99999996247236,9.99999996247236),       // molNStream4 6
    //     mc::Interval(259.16513433898416,259.16513433898416),     // molWStream4 7
    //     mc::Interval(9.994410359108668,9.994410359108668),     // molNStream5 8
    //     mc::Interval(259.01994050480965,259.01994050480965),     // molWStream5 9
    //     mc::Interval(0.005589640891332909,0.005589640891332909),      // molNStream6 10
    //     mc::Interval(0.14519383417451004,0.14519383417451004),     // molWStream6 11
    //     mc::Interval(0.15586818637495115,0.15586818637495115),      // molNStream7 12
    //     mc::Interval(4.039531860542882,4.039531860542882),     // molWStream7 13
    //     mc::Interval(9.83854217273372,9.83854217273372),      // molNStream8 14
    //     mc::Interval(254.98040864426682,254.98040864426682)      // molWStream8 15
    // };
    this->first_stage_IX = {
        mc::Interval(1e-05, 10),      // molNStream1 0
        mc::Interval(1e-05, 500),     // molWStream1 1
        mc::Interval(1e-05, 10),      // molNStream2 2
        mc::Interval(1e-05, 500),     // molWStream2 3
        mc::Interval(1e-05, 10),      // molNStream3 4
        mc::Interval(1e-05, 500),     // molWStream3 5
        mc::Interval(1e-05, 10),      // molNStream4 6
        mc::Interval(1e-05, 500),     // molWStream4 7
        mc::Interval(1e-05, 10),     // molNStream5 8
        mc::Interval(1e-05, 500),     // molWStream5 9
        mc::Interval(1e-05, 10),      // molNStream6 10
        mc::Interval(1e-05, 500),     // molWStream6 11
        mc::Interval(1e-05, 10),      // molNStream7 12
        mc::Interval(1e-05, 500),     // molWStream7 13
        mc::Interval(1e-05, 10),      // molNStream8 14
        mc::Interval(1e-05, 500)      // molWStream8 15
    };

    this->second_stage_IX = {

        mc::Interval(0,1e4),           // currentDensity 16
        mc::Interval(0, 1),           // flowSplit 17
        mc::Interval(1e-2, 10),          // memLength 18
        mc::Interval(1e-2,10),         // memWidth 19
        mc::Interval(0.001, 0.01),     // thicknessConcentrate 20
        mc::Interval(0.001, 0.01),     // thicknessDilute 21
        mc::Interval(1e-3, 100),         // flowOutConcentrateND 22
        mc::Interval(1e-3, 10),         // flowOutDiluteND 23
        mc::Interval(1e-3, 100),          // concOutConcentrateND 24
        mc::Interval(1e-3, 100),          // concOutDiluteND 25
        mc::Interval(0, 5),          // voltCellPair 26
        mc::Interval(1e-6, 1e-2),       // resConcentrate 27
        mc::Interval(1e-6, 1e-2),       // resDilute 28
        mc::Interval(1e-4, 0.25),       // voltNonOhmicCEM 29
        mc::Interval(1e-4, 0.25),       // voltNonOhmicAEM 30
        mc::Interval(0, 1e6),           // capex 31
        mc::Interval(0, 1e6),           // opex1 32
        mc::Interval(0, 1e6),           // opex2 33
        mc::Interval(0, 1e6),            // costTotal 34
        mc::Interval(1e-3, 1e6),            // avgConcDilIntCem 35
        mc::Interval(1e-3, 1e6),            // avgConcDilIntAem 36
    };   
    // this->second_stage_IX = {

    //     mc::Interval(178.74440449245503, 178.74440449245503),     // currentDensity 16
    //     mc::Interval(0.0, 0.0),                                   // flowSplit 17            (raw -1.0e-10, clamped)
    //     mc::Interval(10.0, 10.0),                                 // memLength 18            (raw 10.000000000999169, clamped)
    //     mc::Interval(0.019626763663158974, 0.019626763663158974), // memWidth 19
    //     mc::Interval(0.001, 0.001),                               // thicknessConcentrate 20 (raw 0.0009999999000000574, clamped)
    //     mc::Interval(0.001, 0.001),                               // thicknessDilute 21      (raw 0.0009999999000000063, clamped)
    //     mc::Interval(0.02882349524569438, 0.02882349524569438),   // flowOutConcentrateND 22
    //     mc::Interval(0.7994392490245111, 0.7994392490245111),     // flowOutDiluteND 23
    //     mc::Interval(5.241464713644847, 5.241464713644847),       // concOutConcentrateND 24
    //     mc::Interval(0.07195952395213905, 0.07195952395213905),   // concOutDiluteND 25
    //     mc::Interval(0.42603262990992163, 0.42603262990992163),   // voltCellPair 26
    //     mc::Interval(6.134808984046456e-05, 6.134808984046456e-05), // resConcentrate 27
    //     mc::Interval(0.001040724565921578, 0.001040724565921578), // resDilute 28
    //     mc::Interval(0.07633097808312074, 0.07633097808312074),   // voltNonOhmicCEM 29
    //     mc::Interval(0.07585223733344049, 0.07585223733344049),   // voltNonOhmicAEM 30
    //     mc::Interval(0,1e6),     // capex 31
    //     mc::Interval(0,1e6),     // opex1 32
    //     mc::Interval(0, 1e6),   // opex2 33
    //     mc::Interval(0, 1e6),       // costTotal 34
    //     mc::Interval(1e-3, 1e6),                                  // avgConcDilIntCem 35     (not in Pyomo output — see note)
    //     mc::Interval(1e-3, 1e6),                                  // avgConcDilIntAem 36
    // };       



}

void EDUnits::buildDAG() {
    const double numCells = 500.0;
    const double daysOperation = 7.0 * 365.0;
    const double rg = 8.314;
    const double temp = 273.15 + 25.0;
    const double faraday = 96485.3321;
    const double molweightN = 14.0067;
    const double molweightNh4Cl = 53.491;
    const double thicknessCem = 0.00011;
    const double thicknessAem = 0.000105;
    const double permSelCem = 0.9;
    const double permSelAem = 0.9;
    const double transCem = 0.95;
    const double transAem = 0.95;
    const double transIonConc = 0.491;
    const double transIonDil = 0.509;
    const double saltDiffCem = 1e-12;
    const double saltDiffAem = 1e-12;
    const double saltDiffConc = 1.84e-9;
    const double saltDiffDil = 1.84e-9;
    const double waterPermCem = 7.79e-11 / 3600.0;
    const double waterPermAem = 6.29e-11 / 3600.0;
    const double vantHoffNumber = 2.0;
    const double molweightH20 = 0.018;
    const double molweightH20gMol = 18.01528;
    const double densityH2O = 1000.0;
    const double densityH2OGm3 = 1000.0 * 1000.0;
    const double densityManure = 1000.0;
    const double viscosityManure = 1e-3;
    const double waterTransNumber = 6.0 + 8.0 + 2.0;
    const double osmoticCoeff = 1.0;
    const double resBlank = 2.5e-5;
    const double resCem = 2.5e-4;
    const double resAem = 1.8e-4;
    const double conducConcentrate = 149.6e-4;
    const double conducDilute = 149.6e-4;
    const double costElec = 0.08147;

    const double flowFeed = 192555.0 / (1000.0 * 24.0 * 60.0 * 60.0);
    const double concFeed = 1366.0;
    const double molNFeed = flowFeed * concFeed / molweightN;
    const double molWFeed = flowFeed * (densityH2OGm3 - concFeed) / molweightH20gMol;

    const double scaleFacFlowByUnit[4] = {
        0.0,
        5.388e-6,
        5.388e-6,
        5.388e-6
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        3.967,
        3.967, // please delete e03
        3.967
    };

    for (size_t scenarioIndex = 0; scenarioIndex < this->scenario_names.size(); ++scenarioIndex) {
        auto scenarioName = this->scenario_names[scenarioIndex];
        const int firstStageCount = static_cast<int>(this->first_stage_IX.size());
        const int secondStageCount = static_cast<int>(this->second_stage_IX.size());
        const int nvars = firstStageCount + secondStageCount;

        this->X[scenarioName].resize(nvars);
        for (int i = 0; i < nvars; ++i) {
            this->X[scenarioName][i].set(&this->DAG[scenarioName]);
        }

        auto& vars = this->X[scenarioName];

        // Variable ordering: molNStream1, molWStream1, ..., molNStream8, molWStream8
        mc::FFVar& molNStream1 = vars[0];
        mc::FFVar& molWStream1 = vars[1];
        mc::FFVar& molNStream2 = vars[2];
        mc::FFVar& molWStream2 = vars[3];
        mc::FFVar& molNStream3 = vars[4];
        mc::FFVar& molWStream3 = vars[5];
        mc::FFVar& molNStream4 = vars[6];
        mc::FFVar& molWStream4 = vars[7];
        mc::FFVar& molNStream5 = vars[8];
        mc::FFVar& molWStream5 = vars[9];
        mc::FFVar& molNStream6 = vars[10];
        mc::FFVar& molWStream6 = vars[11];
        mc::FFVar& molNStream7 = vars[12];
        mc::FFVar& molWStream7 = vars[13];
        mc::FFVar& molNStream8 = vars[14];
        mc::FFVar& molWStream8 = vars[15];

        const int secondStageStart = firstStageCount;

        mc::FFVar& currentDensity = vars[secondStageStart + 0];
        mc::FFVar& flowSplit = vars[secondStageStart + 1];
        mc::FFVar& memLength = vars[secondStageStart + 2];
        mc::FFVar& memWidth = vars[secondStageStart + 3];
        mc::FFVar& thicknessConcentrate = vars[secondStageStart + 4];
        mc::FFVar& thicknessDilute = vars[secondStageStart + 5];
        mc::FFVar& flowOutConcentrateNd = vars[secondStageStart + 6];
        mc::FFVar& flowOutDiluteNd = vars[secondStageStart + 7];
        mc::FFVar& concOutConcentrateNd = vars[secondStageStart + 8];
        mc::FFVar& concOutDiluteNd = vars[secondStageStart + 9];
        mc::FFVar& voltCellPair = vars[secondStageStart + 10];
        mc::FFVar& resConcentrate = vars[secondStageStart + 11];
        mc::FFVar& resDilute = vars[secondStageStart + 12];
        mc::FFVar& voltNonOhmicCem = vars[secondStageStart + 13];
        mc::FFVar& voltNonOhmicAem = vars[secondStageStart + 14];
        mc::FFVar& capex = vars[secondStageStart + 15];
        mc::FFVar& opex1 = vars[secondStageStart + 16];
        mc::FFVar& opex2 = vars[secondStageStart + 17];
        mc::FFVar& costTotal = vars[secondStageStart + 18];
        mc::FFVar& avgConcDilIntCem = vars[secondStageStart + 19];
        mc::FFVar& avgConcDilIntAem = vars[secondStageStart + 20];



        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream1*densityH2OGm3 / (molNStream1*molweightN + molWStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream4*densityH2OGm3 / (molNStream4*molweightN + molWStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream5*densityH2OGm3 / (molNStream5*molweightN + molWStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;


        //iSSUES:flowOutConcentrateNd, flowOutDiluteNd,concOutConcentrateNd, concOutDiluteNd

        mc::FFVar reConc = densityManure * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * viscosityManure);
        mc::FFVar reDil = densityManure * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * viscosityManure);

        const double scConc = viscosityManure / (densityManure * saltDiffConc);
        const double scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConc = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDil = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);


        // avgConc*Int* directly, avgConcConc/avgConcDil cancelled out of the correction term
        mc::FFVar avgConcConcIntCem = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc);

        mc::FFVar avgConcConcIntAem = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc);


        // mc::FFVar avgConcDilIntCem = avgConcDil -
        //     (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil);

        // mc::FFVar avgConcDilIntAem = avgConcDil -
        //     (current / (memLength * memWidth)) * (2.0 * thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil); // please delete
        
        std::vector<mc::FFVar> constraints;
        auto addEqualityConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
            constraints.push_back(-expr);
        };
        auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
            constraints.push_back(expr);
        };


        //addEqualityConstraint( avgConcDilIntCem -0.5 * (concOutDiluteNd * scaleFacConc + concInED) +
        //    currentDensity * (2.0 * thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil));

        //addEqualityConstraint( avgConcDilIntAem - 0.5 * (concOutDiluteNd * scaleFacConc + concInED) +
        //    currentDensity * (2.0 * thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil));
        
        
        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (currentDensity / faraday);

        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil)));
        
        
        mc::FFVar fluxIonsTotal = (condFlux + diffFluxAem + diffFluxCem);

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * osmoticCoeff * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * osmoticCoeff * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = (osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux);

        addEqualityConstraint(molNStream1 - molNFeed - molNStream6);
        addEqualityConstraint(molWStream1 - molWFeed - molWStream6);
        addEqualityConstraint(molNStream4 - molNStream3 - molNStream8);
        addEqualityConstraint(molWStream4 - molWStream3 - molWStream8);

        addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1e6) * molWStream2);
        addLessEqualConstraint(30000.0 * (molweightH20gMol / molweightN / 1e6) * molWStream7 - molNStream7);

        addLessEqualConstraint(molNStream1 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream1 * molweightH20gMol);
        addLessEqualConstraint(molNStream2 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream2 * molweightH20gMol);
        addLessEqualConstraint(molNStream3 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream3 * molweightH20gMol);
        addLessEqualConstraint(molNStream4 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream4 * molweightH20gMol);
        addLessEqualConstraint(molNStream5 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream5 * molweightH20gMol);
        addLessEqualConstraint(molNStream6 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream6 * molweightH20gMol);
        addLessEqualConstraint(molNStream7 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream7 * molweightH20gMol);
        addLessEqualConstraint(molNStream8 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream8 * molweightH20gMol);

        addEqualityConstraint(resConcentrate * conducConcentrate * 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * 0.5 * (concOutDiluteNd * scaleFacConc + concInED) - thicknessDilute);


        //addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        //addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);


        addEqualityConstraint(voltNonOhmicCem * (faraday / (permSelCem * rg * temp)) - log(avgConcConcIntCem)+log(avgConcDilIntCem));
        addEqualityConstraint(voltNonOhmicAem * (faraday / (permSelAem * rg * temp)) - log(avgConcConcIntAem) +log(avgConcDilIntAem));


        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));
        

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint(voltCellPair - (voltNonOhmicCem + voltNonOhmicAem) -
            (resConcentrate + resDilute + resCem + resAem) * (currentDensity));


        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * pow(currentDensity, 2.0) + numCells * voltCellPair*currentDensity) * memLength * memWidth * daysOperation);


        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
           densityManure * (1400.0 / reConc) * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate) * ((flowOutConc + flowInConc) / 2.0) +
           densityManure * (104.5 / pow(reDil, 0.37)) * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute) * ((flowOutDil + flowInDil) / 2.0)));

        if (unitIdx == 1) {
            addEqualityConstraint(molNStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 2) {
            addEqualityConstraint(molNStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 3) {
            addEqualityConstraint(molNStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        mc::FFVar objective = this->probability * (capex + opex1 + opex2);

        std::vector<mc::FFVar> functions;

        functions.push_back(objective);
        functions.insert(functions.end(), constraints.begin(), constraints.end());
        this->F[scenarioName] = functions;
    }
}

void EDUnits::buildFullModelDAG() {
    const double numCells = 500.0;
    const double daysOperation = 7.0 * 365.0;
    const double rg = 8.314;
    const double temp = 273.15 + 25.0;
    const double faraday = 96485.3321;
    const double molweightN = 14.0067;
    const double molweightNh4Cl = 53.491;
    const double thicknessCem = 0.00011;
    const double thicknessAem = 0.000105;
    const double permSelCem = 0.9;
    const double permSelAem = 0.9;
    const double transCem = 0.95;
    const double transAem = 0.95;
    const double transIonConc = 0.491;
    const double transIonDil = 0.509;
    const double saltDiffCem = 1e-12;
    const double saltDiffAem = 1e-12;
    const double saltDiffConc = 1.84e-9;
    const double saltDiffDil = 1.84e-9;
    const double waterPermCem = 7.79e-11 / 3600.0;
    const double waterPermAem = 6.29e-11 / 3600.0;
    const double vantHoffNumber = 2.0;
    const double molweightH20 = 0.018;
    const double molweightH20gMol = 18.01528;
    const double densityH2O = 1000.0;
    const double densityH2OGm3 = 1000.0 * 1000.0;
    const double densityManure = 1000.0;
    const double viscosityManure = 1e-3;
    const double waterTransNumber = 6.0 + 8.0 + 2.0;
    const double osmoticCoeff = 1.0;
    const double resBlank = 2.5e-5;
    const double resCem = 2.5e-4;
    const double resAem = 1.8e-4;
    const double conducConcentrate = 149.6e-4;
    const double conducDilute = 149.6e-4;
    const double costElec = 0.08147;

    const double flowFeed = 192555.0 / (1000.0 * 24.0 * 60.0 * 60.0);
    const double concFeed = 1366.0;
    const double molNFeed = flowFeed * concFeed / molweightN;
    const double molWFeed = flowFeed * (densityH2OGm3 - concFeed) / molweightH20gMol;

    const double scaleFacFlowByUnit[4] = {
        0.0,
        5.388e-6,
        5.388e-6,
        5.388e-6
    };
    const double scaleFacConcByUnit[4] = {
        0.0,
        3.967e-5,
        3.967e-5,
        3.967e-5
    };

    const int firstStageCount = static_cast<int>(this->first_stage_IX.size());
    // second_stage_IX is already the full combined size after convertToCentralizedModel()
    const int totalSecondStageCount = static_cast<int>(this->second_stage_IX.size());
    const int secondStageCount = totalSecondStageCount / static_cast<int>(this->scenario_names.size());

    this->X[ScenarioNames::SCENARIO1].resize(firstStageCount + totalSecondStageCount);
    for (int i = 0; i < firstStageCount + totalSecondStageCount; ++i) {
        this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);
    }

    auto& vars = this->X[ScenarioNames::SCENARIO1];

    mc::FFVar& molNStream1 = vars[0];
    mc::FFVar& molWStream1 = vars[1];
    mc::FFVar& molNStream2 = vars[2];
    mc::FFVar& molWStream2 = vars[3];
    mc::FFVar& molNStream3 = vars[4];
    mc::FFVar& molWStream3 = vars[5];
    mc::FFVar& molNStream4 = vars[6];
    mc::FFVar& molWStream4 = vars[7];
    mc::FFVar& molNStream5 = vars[8];
    mc::FFVar& molWStream5 = vars[9];
    mc::FFVar& molNStream6 = vars[10];
    mc::FFVar& molWStream6 = vars[11];
    mc::FFVar& molNStream7 = vars[12];
    mc::FFVar& molWStream7 = vars[13];
    mc::FFVar& molNStream8 = vars[14];
    mc::FFVar& molWStream8 = vars[15];


    std::vector<mc::FFVar> constraints;


    auto addEqualityConstraint = [&](const mc::FFVar& expr) {
        constraints.push_back(expr);
        constraints.push_back(-expr);
    };
    auto addLessEqualConstraint = [&](const mc::FFVar& expr) {
        constraints.push_back(expr);
    };

    addEqualityConstraint(molNStream1 - molNFeed - molNStream6);
    addEqualityConstraint(molWStream1 - molWFeed - molWStream6);
    addEqualityConstraint(molNStream4 - molNStream3 - molNStream8);
    addEqualityConstraint(molWStream4 - molWStream3 - molWStream8);

    addLessEqualConstraint(molNStream2 - 400.0 * (molweightH20gMol / molweightN / 1e6) * molWStream2);
    addLessEqualConstraint(30000.0 * (molweightH20gMol / molweightN / 1e6) * molWStream7 - molNStream7);

    addLessEqualConstraint(molNStream1 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream1 * molweightH20gMol);
    addLessEqualConstraint(molNStream2 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream2 * molweightH20gMol);
    addLessEqualConstraint(molNStream3 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream3 * molweightH20gMol);
    addLessEqualConstraint(molNStream4 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream4 * molweightH20gMol);
    addLessEqualConstraint(molNStream5 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream5 * molweightH20gMol);
    addLessEqualConstraint(molNStream6 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream6 * molweightH20gMol);
    addLessEqualConstraint(molNStream7 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream7 * molweightH20gMol);
    addLessEqualConstraint(molNStream8 * (densityH2OGm3 - 40000*molweightN) - 40000 * molWStream8 * molweightH20gMol);

    mc::FFVar objective = 0.0;

    for (size_t scenarioIndex = 0; scenarioIndex < this->scenario_names.size(); ++scenarioIndex) {
        const int unitIdx = static_cast<int>(scenarioIndex) + 1;
        const int secondStageStart = firstStageCount + static_cast<int>(scenarioIndex) * secondStageCount;

        mc::FFVar& currentDensity = vars[secondStageStart + 0];
        mc::FFVar& flowSplit = vars[secondStageStart + 1];
        mc::FFVar& memLength = vars[secondStageStart + 2];
        mc::FFVar& memWidth = vars[secondStageStart + 3];
        mc::FFVar& thicknessConcentrate = vars[secondStageStart + 4];
        mc::FFVar& thicknessDilute = vars[secondStageStart + 5];
        mc::FFVar& flowOutConcentrateNd = vars[secondStageStart + 6];
        mc::FFVar& flowOutDiluteNd = vars[secondStageStart + 7];
        mc::FFVar& concOutConcentrateNd = vars[secondStageStart + 8];
        mc::FFVar& concOutDiluteNd = vars[secondStageStart + 9];
        mc::FFVar& voltCellPair = vars[secondStageStart + 10];
        mc::FFVar& resConcentrate = vars[secondStageStart + 11];
        mc::FFVar& resDilute = vars[secondStageStart + 12];
        mc::FFVar& voltNonOhmicCem = vars[secondStageStart + 13];
        mc::FFVar& voltNonOhmicAem = vars[secondStageStart + 14];
        mc::FFVar& capex = vars[secondStageStart + 15];
        mc::FFVar& opex1 = vars[secondStageStart + 16];
        mc::FFVar& opex2 = vars[secondStageStart + 17];
        mc::FFVar& costTotal = vars[secondStageStart + 18];
        mc::FFVar& avgConcDilIntCem = vars[secondStageStart + 19];
        mc::FFVar& avgConcDilIntAem = vars[secondStageStart + 20];



        const double scaleFacFlow = scaleFacFlowByUnit[unitIdx];
        const double scaleFacConc = scaleFacConcByUnit[unitIdx];

        mc::FFVar flowInED, concInED;
        if (unitIdx == 1) {
            flowInED = (molNStream1 * molweightN + molWStream1 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream1*densityH2OGm3 / (molNStream1*molweightN + molWStream1 * molweightH20gMol);
        } else if (unitIdx == 2) {
            flowInED = (molNStream4 * molweightN + molWStream4 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream4*densityH2OGm3 / (molNStream4*molweightN + molWStream4 * molweightH20gMol);
        } else {
            flowInED = (molNStream5 * molweightN + molWStream5 * molweightH20gMol) / densityH2OGm3;
            concInED = molNStream5*densityH2OGm3 / (molNStream5*molweightN + molWStream5 * molweightH20gMol);
        }

        mc::FFVar flowInConc = flowSplit * flowInED / numCells;
        mc::FFVar flowInDil = (1.0 - flowSplit) * flowInED / numCells;

        mc::FFVar reConc = densityManure * (flowOutConcentrateNd * scaleFacFlow + flowInConc) / (memWidth * viscosityManure);
        mc::FFVar reDil = densityManure * (flowOutDiluteNd * scaleFacFlow + flowInDil) / (memWidth * viscosityManure);

        const double scConc = viscosityManure / (densityManure * saltDiffConc);
        const double scDil = viscosityManure / (densityManure * saltDiffDil);

        mc::FFVar shConc = 0.29 * pow(reConc, 0.5) * pow(scConc, 0.33);
        mc::FFVar shDil = 0.29 * pow(reDil, 0.5) * pow(scDil, 0.33);

        // avgConc*Int* directly, avgConcConc/avgConcDil cancelled out of the correction term
        mc::FFVar avgConcConcIntCem = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc);

        mc::FFVar avgConcConcIntAem = 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc);


        addEqualityConstraint( avgConcDilIntCem -0.5 * (concOutDiluteNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil));

        addEqualityConstraint( avgConcDilIntAem - 0.5 * (concOutDiluteNd * scaleFacConc + concInED) +
            currentDensity * (2.0 * thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil));


        mc::FFVar condFlux = (transCem - (1.0 - transAem)) * (currentDensity / faraday);

        mc::FFVar diffFluxAem = -(saltDiffAem / thicknessAem) * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar diffFluxCem = -(saltDiffCem / thicknessCem) * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar fluxIonsTotal = condFlux + diffFluxAem + diffFluxCem;

        mc::FFVar osmWaterFluxAem = waterPermAem * vantHoffNumber * rg * temp * osmoticCoeff * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transAem - transIonDil)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transAem - transIonDil)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar osmWaterFluxCem = waterPermCem * vantHoffNumber * rg * temp * osmoticCoeff * (0.5 * (concOutConcentrateNd * scaleFacConc-concOutDiluteNd * scaleFacConc) +
            2*currentDensity * ((thicknessConcentrate * (transCem - transIonConc)) / (shConc * faraday * saltDiffConc)
                            +(thicknessDilute * (transCem - transIonConc)) / (shDil * faraday * saltDiffDil)));

        mc::FFVar eosmWaterFlux = waterTransNumber * fluxIonsTotal * molweightH20 / densityH2O;
        mc::FFVar fluxWaterTotal = osmWaterFluxAem + osmWaterFluxCem + eosmWaterFlux;

        addEqualityConstraint(resConcentrate * conducConcentrate * 0.5 * (concOutConcentrateNd * scaleFacConc + concInED) - thicknessConcentrate);
        addEqualityConstraint(resDilute * conducDilute * 0.5 * (concOutDiluteNd * scaleFacConc + concInED) - thicknessDilute);


        //addEqualityConstraint(exp(voltNonOhmicCem * faraday / (permSelCem * rg * temp)) * avgConcDilIntCem - avgConcConcIntCem);
        //addEqualityConstraint(exp(voltNonOhmicAem * faraday / (permSelAem * rg * temp)) * avgConcDilIntAem - avgConcConcIntAem);
        addEqualityConstraint(voltNonOhmicCem * (faraday / (permSelCem * rg * temp)) - log(avgConcConcIntCem)+log(avgConcDilIntCem));

        addEqualityConstraint(voltNonOhmicAem * (faraday / (permSelAem * rg * temp)) - log(avgConcConcIntAem) +log(avgConcDilIntAem));

        addEqualityConstraint(flowOutConcentrateNd - (flowInConc / scaleFacFlow) -
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) -
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint(flowOutDiluteNd - (flowInDil / scaleFacFlow) +
            (fluxWaterTotal * memLength * memWidth / scaleFacFlow) +
            ((molweightNh4Cl * 0.001 / densityH2O) * fluxIonsTotal * memLength * memWidth / scaleFacFlow));

        addEqualityConstraint((flowOutConcentrateNd * concOutConcentrateNd) -
            (flowInConc * concInED) / (scaleFacFlow * scaleFacConc) -
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));

        addEqualityConstraint((flowOutDiluteNd * concOutDiluteNd) -
            (flowInDil * concInED) / (scaleFacFlow * scaleFacConc) +
            fluxIonsTotal * memLength * memWidth / (scaleFacFlow * scaleFacConc));


        addEqualityConstraint(voltCellPair - (voltNonOhmicCem + voltNonOhmicAem) -
            (resConcentrate + resDilute + resCem + resAem) * (currentDensity));

        addEqualityConstraint(capex - (6800.0 * memLength * memWidth + 2.0 * 100.0 * (numCells - 2.0) * memLength * memWidth));

        addEqualityConstraint(opex1 - costElec * (24.0 * 1e-3) *
            (resBlank * pow(currentDensity, 2.0) + numCells * voltCellPair*currentDensity) * memLength * memWidth * daysOperation);

        mc::FFVar flowOutConc = flowOutConcentrateNd * scaleFacFlow;
        mc::FFVar flowOutDil = flowOutDiluteNd * scaleFacFlow;
        addEqualityConstraint(opex2 - costElec * (24.0 * 1e-3) * numCells * daysOperation * (
            densityManure * (1400.0 / reConc) * memLength * pow(flowOutConc / memWidth, 2.0) / (4.0 * thicknessConcentrate) * ((flowOutConc + flowInConc) / 2.0) +
            densityManure * (104.5 / pow(reDil, 0.37)) * memLength * pow(flowOutDil / memWidth, 2.0) / (4.0 * thicknessDilute) * ((flowOutDil + flowInDil) / 2.0)));

        if (unitIdx == 1) {
            addEqualityConstraint(molNStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream2 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream3 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 2) {
            addEqualityConstraint(molNStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream5 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream6 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        if (unitIdx == 3) {
            addEqualityConstraint(molNStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * concOutConcentrateNd * scaleFacConc);
            addEqualityConstraint(molWStream7 - flowOutConcentrateNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutConcentrateNd * scaleFacConc * molweightN) / molweightH20gMol);
            addEqualityConstraint(molNStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * concOutDiluteNd * scaleFacConc);
            addEqualityConstraint(molWStream8 - flowOutDiluteNd * (scaleFacFlow * numCells) * (densityH2OGm3 - concOutDiluteNd * scaleFacConc * molweightN) / molweightH20gMol);
        }

        objective += this->probability * (capex + opex1 + opex2);
    }

    std::vector<mc::FFVar> functions;

    functions.push_back(objective);
    functions.insert(functions.end(), constraints.begin(), constraints.end());
    this->F[ScenarioNames::SCENARIO1] = functions;
    this->full_model_built = true;
}

Ipopt::SmartPtr<STModel> EDUnits::clone() {
    Ipopt::SmartPtr<EDUnits> copy = new EDUnits();

    copy->scenario_name = this->scenario_name;
    copy->first_stage_IX = this->first_stage_IX;
    copy->second_stage_IX = this->second_stage_IX;
    copy->scenario_names = this->scenario_names;
    copy->probability = this->probability;
    copy->clearDAG();
    if (this->full_model_built) {
        copy->buildFullModelDAG();
    } else {
        copy->buildDAG();
    }
    return copy;
}