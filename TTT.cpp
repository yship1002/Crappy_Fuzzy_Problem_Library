#include "TTT.h"

TTT::TTT(BranchingStrategy branching_strategy):STModel() {

    this->branching_strategy = branching_strategy;
    this->scenario_names = { ScenarioNames::SCENARIO1,
        ScenarioNames::SCENARIO2,ScenarioNames::SCENARIO3,ScenarioNames::SCENARIO4,
        ScenarioNames::SCENARIO5,ScenarioNames::SCENARIO6,ScenarioNames::SCENARIO7,ScenarioNames::SCENARIO8,
        ScenarioNames::SCENARIO9,ScenarioNames::SCENARIO10,ScenarioNames::SCENARIO11,ScenarioNames::SCENARIO12,
        ScenarioNames::SCENARIO13,ScenarioNames::SCENARIO14,ScenarioNames::SCENARIO15,ScenarioNames::SCENARIO16,
        ScenarioNames::SCENARIO17,ScenarioNames::SCENARIO18,ScenarioNames::SCENARIO19,ScenarioNames::SCENARIO20,
        ScenarioNames::SCENARIO21,ScenarioNames::SCENARIO22,ScenarioNames::SCENARIO23,ScenarioNames::SCENARIO24,
        ScenarioNames::SCENARIO25,ScenarioNames::SCENARIO26,ScenarioNames::SCENARIO27,ScenarioNames::SCENARIO28,
        ScenarioNames::SCENARIO29,ScenarioNames::SCENARIO30,ScenarioNames::SCENARIO31,ScenarioNames::SCENARIO32,
        ScenarioNames::SCENARIO33,ScenarioNames::SCENARIO34,ScenarioNames::SCENARIO35,ScenarioNames::SCENARIO36,
        ScenarioNames::SCENARIO37,ScenarioNames::SCENARIO38,ScenarioNames::SCENARIO39,ScenarioNames::SCENARIO40,
        ScenarioNames::SCENARIO41,ScenarioNames::SCENARIO42,ScenarioNames::SCENARIO43,ScenarioNames::SCENARIO44,
        ScenarioNames::SCENARIO45,ScenarioNames::SCENARIO46,ScenarioNames::SCENARIO47,ScenarioNames::SCENARIO48,
        ScenarioNames::SCENARIO49,ScenarioNames::SCENARIO50,ScenarioNames::SCENARIO51,ScenarioNames::SCENARIO52,
        ScenarioNames::SCENARIO53,ScenarioNames::SCENARIO54,ScenarioNames::SCENARIO55,ScenarioNames::SCENARIO56,
        ScenarioNames::SCENARIO57,ScenarioNames::SCENARIO58,ScenarioNames::SCENARIO59,ScenarioNames::SCENARIO60,
        ScenarioNames::SCENARIO61,ScenarioNames::SCENARIO62,ScenarioNames::SCENARIO63,ScenarioNames::SCENARIO64,
        ScenarioNames::SCENARIO65,ScenarioNames::SCENARIO66,ScenarioNames::SCENARIO67,ScenarioNames::SCENARIO68,
        ScenarioNames::SCENARIO69,ScenarioNames::SCENARIO70,ScenarioNames::SCENARIO71,ScenarioNames::SCENARIO72,
        ScenarioNames::SCENARIO73,ScenarioNames::SCENARIO74,ScenarioNames::SCENARIO75,ScenarioNames::SCENARIO76,
        ScenarioNames::SCENARIO77,ScenarioNames::SCENARIO78,ScenarioNames::SCENARIO79,ScenarioNames::SCENARIO80,
        ScenarioNames::SCENARIO81,ScenarioNames::SCENARIO82,ScenarioNames::SCENARIO83,ScenarioNames::SCENARIO84,
        ScenarioNames::SCENARIO85,ScenarioNames::SCENARIO86,ScenarioNames::SCENARIO87,ScenarioNames::SCENARIO88,
        ScenarioNames::SCENARIO89,ScenarioNames::SCENARIO90,ScenarioNames::SCENARIO91,ScenarioNames::SCENARIO92,
        ScenarioNames::SCENARIO93,ScenarioNames::SCENARIO94,ScenarioNames::SCENARIO95,ScenarioNames::SCENARIO96,
        ScenarioNames::SCENARIO97,ScenarioNames::SCENARIO98,ScenarioNames::SCENARIO99,ScenarioNames::SCENARIO100

    };
    this->scenario_name = ScenarioNames::SCENARIO1; //default
    this->probability = 1; // equal probability for each scenario

    // feed-type driver -> e5's bound directly, same values as before

    std::vector<double> fix_1;
    std::vector<double> fix_2;

    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            fix_1.push_back(-5 + ((-1) - (-5)) * static_cast<double>(i) / 9);
            fix_2.push_back(2 + ((6) - (2)) * static_cast<double>(j) / 9);
        }

    }
    for (auto& scenario_name : this->scenario_names) {
        this->perturb_coeffs[scenario_name] = std::vector<double>{fix_1[static_cast<int>(scenario_name)], fix_2[static_cast<int>(scenario_name)]};
    }



    // NOTE: roles now match the Python const_model() (4 first-stage
    // variables x0-x3, bounds [0,1], shared across scenarios / no
    // scenario subscript; 2 second-stage variables x4,x5 per scenario,
    // bounds [1e-5,16]).
    this->first_stage_IX = {
        mc::Interval(-5,-1), //x0
        mc::Interval(2,6) //x1
    };

    // Per-scenario second-stage bounds (x4, x5). This holds ONE
    // scenario's worth of second-stage variables, not the total across
    // all scenarios -- buildFullModelDAG() must multiply this size by
    // scenario_names.size() to get the full-model variable count.
    this->second_stage_IX = {

    };

};
void TTT::buildDAG() {
    for (auto& scenario_name : this->scenario_names) {
        int n_first_stage_vars = this->first_stage_IX.size();

        // Loop over each scenario to build subproblem

        const int nvars = n_first_stage_vars + this->second_stage_IX.size();

        this->X[scenario_name].resize(nvars);


        for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);

        // scenario perturbation and coefficients


        mc::FFVar c1=this->X[scenario_name][0]-this->perturb_coeffs[scenario_name][0];
        mc::FFVar nc1=-c1;
        mc::FFVar c2=this->X[scenario_name][1]-this->perturb_coeffs[scenario_name][1];
        mc::FFVar nc2=-c2;

        // objective: prob*(-1e6*x3 + recourse_cost[s]*(x4+x5)), matching
        // Python's _obj exactly.
        mc::FFVar objective = -this->X[scenario_name][0]/this->X[scenario_name][1];
        this->F[scenario_name]={objective,c1,nc1,c2,nc2};
    }
}
void TTT::buildFullModelDAG(){

}
Ipopt::SmartPtr<STModel> TTT::clone(){
    Ipopt::SmartPtr<TTT> p = new TTT();

    p->scenario_name=this->scenario_name;
    p->first_stage_IX=this->first_stage_IX;
    p->second_stage_IX=this->second_stage_IX;
    p->perturb=this->perturb;
    p->perturb_coeffs=this->perturb_coeffs;

    p->scenario_names=this->scenario_names;
    p->probability=this->probability;
    p->first_stage_map=this->first_stage_map;
    p->second_stage_map=this->second_stage_map;
    p->clearDAG(); // clear the DAG of the cloned model
    if (this->full_model_built) {
        p->buildFullModelDAG();
    } else {
        p->buildDAG();
    }
    return p;
}