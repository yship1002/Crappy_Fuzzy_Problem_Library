#include "Ex722.h"

Ex722Model::Ex722Model(BranchingStrategy branching_strategy, int num_scenarios):STModel() {

    this->branching_strategy = branching_strategy;
    this->scenario_names = { ScenarioNames::SCENARIO1, ScenarioNames::SCENARIO2,ScenarioNames::SCENARIO3,
        ScenarioNames::SCENARIO4, ScenarioNames::SCENARIO5
    };
    this->scenario_name = ScenarioNames::SCENARIO1; //default
    this->probability = 0.2; // equal probability for each scenario

    // ===== ORIGINAL Ex722 scenarios (near-identical; UBD ~ -378487.78) =====
    // To REVERT: uncomment this block and comment out the NEW block below.
    // // feed-type driver -> e5's bound directly, same values as before
    // this->perturb = {
    //     {ScenarioNames::SCENARIO1, 10.0},
    //     {ScenarioNames::SCENARIO2, 20.0},
    //     {ScenarioNames::SCENARIO3, 30.0},
    //     {ScenarioNames::SCENARIO4, 40.0},
    //     {ScenarioNames::SCENARIO5, 15}
    // };
    //
    // // ------------------------------------------------------------------
    // // Scenario data matching the Python const_model(): explicit lists,
    // // one value per scenario, for the underlying uncertain quantities.
    // // temp_factor scales the bilinear coupling coefficients (a1-a4);
    // // price_shock scales the second-stage recourse cost. Every
    // // scenario-dependent coefficient below is derived from these lists,
    // // in the same order as this->scenario_names.
    // // ------------------------------------------------------------------
    // std::vector<double> temp_factor = {
    //     0.85, 0.90, 0.95, 1.00, 1.05
    // };
    // std::vector<double> price_shock = {
    //     1.00, 1.05, 0.95, 1.10, 0.90
    // };
    //
    // for (size_t i = 0; i < this->scenario_names.size(); ++i) {
    //     ScenarioNames sn = this->scenario_names[i];
    //     double tf = temp_factor[i];
    //
    //     // {a1, a2, a3, a4} for this scenario, matching Python's
    //     // a1_vals/a2_vals/a3_vals/a4_vals
    //     this->perturb_coeffs[sn] = {
    //         0.09755988 * tf,             // a1
    //         0.0965842812 * std::sqrt(tf),// a2
    //         0.0391908 * tf,              // a3
    //         0.03527172 * std::sqrt(tf)   // a4
    //     };
    //
    //     this->recourse_cost[sn] = 1000.0 * price_shock[i];
    // }

    // ===== NEW diverse 5-scenario problem (UBD -92551.4) =====
    // BEGIN NEW BLOCK
    // ------------------------------------------------------------------
    // Five deliberately DIFFERENT scenarios (were near-copies before).
    //
    // Structure: the equalities e1-e4 only keep a common first-stage
    // manifold across scenarios if the ratios a3/a1 and a4/a2 are the
    // same in every scenario (otherwise no x0..x3 is feasible for all
    // scenarios). So each scenario scales the pair (a1,a3) by k1 and
    // the pair (a2,a4) by k2 over orders of magnitude; this moves the
    // recourse values x4 = (1-x0)/(a1*x1), x5 = (x0-x1)/(a2*x2) and
    // hence the e5 (sqrt budget) / x4,x5-bound cuts on the first stage.
    // perturb[s] is the e5 budget and recourse_cost[s] the recourse price.
    //
    //   scen  k1    k2    budget  recourse   character
    //   S1    3.5   0.4   12.5    45000      costly recourse, mid x3
    //   S2    0.4   5.0   5.5     40000      tight budget + costly: wants x3~0.1
    //   S3    0.8   0.4   15.5    25         cheap recourse: wants x3 max (~0.40)
    //   S4    1.0   15.0  9.0     100        cheap, x5 tiny: wants x3 max
    //   S5    0.8   2.5   7.5     60000      costly recourse: wants x3~0.1
    //
    // Individually the scenarios prefer clearly different first-stage
    // points (S3,S4: x3~0.40; S1: x3~0.33; S2,S5: x3~0.10) and their
    // feasible first-stage sets cover 0.3%-3.5% of the grid, yet the
    // intersection is non-empty (checked numerically).
    // ------------------------------------------------------------------
    this->perturb = {
        {ScenarioNames::SCENARIO1, 12.5},
        {ScenarioNames::SCENARIO2, 5.5},
        {ScenarioNames::SCENARIO3, 15.5},
        {ScenarioNames::SCENARIO4, 9.0},
        {ScenarioNames::SCENARIO5, 7.5}
    };
    std::vector<double> k1 = {3.5, 0.4, 0.8, 1.0, 0.8};   // scales a1,a3
    std::vector<double> k2 = {0.4, 5.0, 0.4, 15.0, 2.5};  // scales a2,a4
    std::vector<double> rcs = {45000.0, 40000.0, 25.0, 100.0, 60000.0};

    for (size_t i = 0; i < this->scenario_names.size(); ++i) {
        ScenarioNames sn = this->scenario_names[i];
        // {a1, a2, a3, a4}; base values from the Python const_model()
        this->perturb_coeffs[sn] = {
            0.09755988 * k1[i],     // a1
            0.0965842812 * k2[i],   // a2
            0.0391908 * k1[i],      // a3 (a3/a1 fixed)
            0.03527172 * k2[i]      // a4 (a4/a2 fixed)
        };
        this->recourse_cost[sn] = rcs[i];
    }
    // END NEW BLOCK

    // ===== extra scenarios 6-20 (setScenarioCount() below keeps the first num_scenarios) =====
    // 15 extra scenarios, all OUTSIDE/cross-paired against the original 5
    // (original envelope: k1 [0.4,3.5], k2 [0.4,15], budget [5.5,15.5],
    // recourse [25,60000]; original pattern: costly recourse <-> tight
    // budget, cheap <-> loose). New ones use k1 [0.15,5], k2 [0.15,20],
    // budget [4.5,20], recourse [3,95000], and include costly+loose and
    // cheap+tight pairings. a3/a1 and a4/a2 ratios stay fixed (see above).
    {
        const ScenarioNames new_names[15] = {
            ScenarioNames::SCENARIO6,  ScenarioNames::SCENARIO7,  ScenarioNames::SCENARIO8,
            ScenarioNames::SCENARIO9,  ScenarioNames::SCENARIO10, ScenarioNames::SCENARIO11,
            ScenarioNames::SCENARIO12, ScenarioNames::SCENARIO13, ScenarioNames::SCENARIO14,
            ScenarioNames::SCENARIO15, ScenarioNames::SCENARIO16, ScenarioNames::SCENARIO17,
            ScenarioNames::SCENARIO18, ScenarioNames::SCENARIO19, ScenarioNames::SCENARIO20};
        const double new_k1[15]  = {5.0, 0.15, 2.0, 0.2, 4.0, 0.3, 1.5, 2.5, 0.15, 2.5, 0.6, 3.0, 1.2, 0.25, 4.5};
        const double new_k2[15]  = {0.15, 20.0, 8.0, 0.2, 20.0, 1.0, 0.15, 0.3, 3.0, 0.6, 12.0, 0.5, 0.3, 18.0, 1.5};
        const double new_bud[15] = {18.0, 4.8, 20.0, 6.0, 12.0, 17.0, 5.0, 8.0, 19.0, 4.5, 14.0, 6.5, 11.0, 10.0, 16.0};
        const double new_rcs[15] = {80000.0, 10.0, 70000.0, 5.0, 90000.0, 30000.0, 8.0, 50.0, 85000.0, 15.0, 75000.0, 12.0, 95000.0, 3.0, 40000.0};
        for (int j = 0; j < 15; ++j) {
            ScenarioNames sn = new_names[j];
            this->scenario_names.push_back(sn);
            this->perturb[sn] = new_bud[j];
            this->perturb_coeffs[sn] = {
                0.09755988 * new_k1[j],
                0.0965842812 * new_k2[j],
                0.0391908 * new_k1[j],
                0.03527172 * new_k2[j]
            };
            this->recourse_cost[sn] = new_rcs[j];
        }
    }
    // ===== end extra scenarios =====
    // Equal weights; with the block above commented out this is 0.2 again.
    // (original literal: this->probability = 0.2;)
    this->setScenarioCount(num_scenarios); // keep first num_scenarios of the 20 defined above
    trimScenarioMap(this->recourse_cost, num_scenarios);
    this->probability = 1.0 / this->scenario_names.size();

    // NOTE: roles now match the Python const_model() (4 first-stage
    // variables x0-x3, bounds [0,1], shared across scenarios / no
    // scenario subscript; 2 second-stage variables x4,x5 per scenario,
    // bounds [1e-5,16]).
    this->first_stage_IX = {
        mc::Interval(0, 1), //x0
        mc::Interval(0, 1), //x1
        mc::Interval(0, 1), //x2
        mc::Interval(0, 1)  //x3
    };

    // Per-scenario second-stage bounds (x4, x5). This holds ONE
    // scenario's worth of second-stage variables, not the total across
    // all scenarios -- buildFullModelDAG() must multiply this size by
    // scenario_names.size() to get the full-model variable count.
    this->second_stage_IX = {
        mc::Interval(1e-5,16), //x4
        mc::Interval(1e-5,16)  //x5
    };

};
void Ex722Model::buildDAG() {
    // Only build the scenario currently being asked for (this->scenario_name)
    {
        const ScenarioNames scenario_name = this->scenario_name;
        int n_first_stage_vars = this->first_stage_IX.size();

        // Loop over each scenario to build subproblem

        const int nvars = n_first_stage_vars + this->second_stage_IX.size();

        this->X[scenario_name].resize(nvars);


        for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);

        // scenario perturbation and coefficients
        double p = this->perturb[scenario_name];
        const auto& a = this->perturb_coeffs[scenario_name];
        double a1 = a[0], a2 = a[1], a3 = a[2], a4 = a[3];
        double rc = this->recourse_cost[scenario_name];

        mc::FFVar c1,c2,c3,c4,c5;
        mc::FFVar nc1,nc2,nc3,nc4;

        // Variable layout for this scenario: X[0..3] = x0..x3 (first
        // stage), X[4],X[5] = x4,x5 (second stage). Matches Python's
        // e1-e5 exactly (index-for-index, not just up to relabeling).

        //.  a1*x1*x4 + x0 == 1;
        c1=a1*this->X[scenario_name][1]*this->X[scenario_name][4]+this->X[scenario_name][0]-1;
        nc1=-c1;

        //.  a2*x2*x5 + x1 - x0 == 0;
        c2=a2*this->X[scenario_name][2]*this->X[scenario_name][5]+this->X[scenario_name][1]-this->X[scenario_name][0];
        nc2=-c2;

        //.  a3*x2*x4 + x2 + x0 == 1;
        c3=a3*this->X[scenario_name][2]*this->X[scenario_name][4]+this->X[scenario_name][2]+this->X[scenario_name][0]-1;
        nc3=-c3;

        // a4*x3*x5 + x3 - x0 + x1 - x2 == 0;
        c4=a4*this->X[scenario_name][3]*this->X[scenario_name][5]+this->X[scenario_name][3]-this->X[scenario_name][0]+this->X[scenario_name][1]-this->X[scenario_name][2];
        nc4=-c4;

        // x4**0.5 + x5**0.5 =L= perturb[s];
        c5=pow(this->X[scenario_name][4], 0.5)+pow(this->X[scenario_name][5], 0.5)-p;

        // objective: prob*(-1e6*x3 + recourse_cost[s]*(x4+x5)), matching
        // Python's _obj exactly.
        mc::FFVar objective = -1000000*this->probability*this->X[scenario_name][3]
                             + this->probability*rc*(this->X[scenario_name][4]+this->X[scenario_name][5]);
        this->F[scenario_name]={objective,c1,c2,c3,c4,c5,nc1,nc2,nc3,nc4};
    }
}
void Ex722Model::buildFullModelDAG(){
    // for full model solve we will stay in scenario 1
    int n_first_stage_vars = this->first_stage_IX.size();

    // FIX: second_stage_IX holds ONE scenario's second-stage vars
    // (x4, x5 -> size 2). It must NOT be divided by scenario_names.size();
    // that previously gave n_second_stage_vars = 2/5 = 0 (integer
    // division), which undersized the allocation and made every
    // scenario silently alias the same two "second stage" slots.
    int n_second_stage_vars = this->second_stage_IX.size()/this->scenario_names.size(); // per-scenario count (2: x4, x5)
    int nvars = n_first_stage_vars + n_second_stage_vars * this->scenario_names.size(); // 4 + 2*5 = 14

    this->X[ScenarioNames::SCENARIO1].resize(nvars);

    for (int i = 0; i < n_first_stage_vars; ++i) this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);
    for (int s_idx=0; s_idx<this->scenario_names.size(); ++s_idx){
        int second_stage_start_idx = n_first_stage_vars + s_idx * n_second_stage_vars;
        for (int i = 0; i < n_second_stage_vars; ++i){
            this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+i].set(&this->DAG[ScenarioNames::SCENARIO1]);
        }
    }
    mc::FFVar objective=0;
    for (int s_idx=0; s_idx<this->scenario_names.size(); ++s_idx){
        int second_stage_start_idx = n_first_stage_vars + s_idx * n_second_stage_vars;

        // scenario perturbation and coefficients
        ScenarioNames sn = this->scenario_names[s_idx];
        double p = this->perturb[sn];
        const auto& a = this->perturb_coeffs[sn];
        double a1 = a[0], a2 = a[1], a3 = a[2], a4 = a[3];
        double rc = this->recourse_cost[sn];

        mc::FFVar c1,c2,c3,c4,c5;
        mc::FFVar nc1,nc2,nc3,nc4;

        // X[0..3] = x0..x3 (first stage, shared across scenarios, no
        // scenario subscript - matches Python where x0-x3 are plain
        // Vars, not indexed by pm.S). X[second_stage_start_idx],
        // X[second_stage_start_idx+1] = this scenario's x4[s], x5[s].

        //.  a1*x1*x4[s] + x0 == 1;
        c1=a1*this->X[ScenarioNames::SCENARIO1][1]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][0]-1;
        nc1=-c1;

        //.  a2*x2*x5[s] + x1 - x0 == 0;
        c2=a2*this->X[ScenarioNames::SCENARIO1][2]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]+this->X[ScenarioNames::SCENARIO1][1]-this->X[ScenarioNames::SCENARIO1][0];
        nc2=-c2;

        //.  a3*x2*x4[s] + x2 + x0 == 1;
        c3=a3*this->X[ScenarioNames::SCENARIO1][2]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][2]+this->X[ScenarioNames::SCENARIO1][0]-1;
        nc3=-c3;

        // a4*x3*x5[s] + x3 - x0 + x1 - x2 == 0;
        c4=a4*this->X[ScenarioNames::SCENARIO1][3]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]+this->X[ScenarioNames::SCENARIO1][3]-this->X[ScenarioNames::SCENARIO1][0]+this->X[ScenarioNames::SCENARIO1][1]-this->X[ScenarioNames::SCENARIO1][2];
        nc4=-c4;

        // x4[s]**0.5 + x5[s]**0.5 =L= perturb[s];
        c5=pow(this->X[ScenarioNames::SCENARIO1][second_stage_start_idx], 0.5)+pow(this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1], 0.5)-p;

        objective += -1000000*this->probability*this->X[ScenarioNames::SCENARIO1][3]
                   + this->probability*rc*(this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]);

        std::vector<mc::FFVar> scenario_constraints = {c1,c2,c3,c4,nc1,nc2,nc3,nc4,c5};
        this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].end(), scenario_constraints.begin(), scenario_constraints.end());
    }
    this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].begin(), objective);
    this->full_model_built = true;
}
Ipopt::SmartPtr<STModel> Ex722Model::clone(){
    Ipopt::SmartPtr<Ex722Model> p = new Ex722Model();

    p->scenario_name=this->scenario_name;
    p->first_stage_IX=this->first_stage_IX;
    p->second_stage_IX=this->second_stage_IX;
    p->perturb=this->perturb;
    p->perturb_coeffs=this->perturb_coeffs;
    p->recourse_cost=this->recourse_cost;
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