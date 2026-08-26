// #include "Ex722.h"

// Ex722Model::Ex722Model(BranchingStrategy branching_strategy):STModel() {

//     this->branching_strategy = branching_strategy;
//     this->scenario_names = { ScenarioNames::SCENARIO1, ScenarioNames::SCENARIO2,ScenarioNames::SCENARIO3,
//         ScenarioNames::SCENARIO4, ScenarioNames::SCENARIO5,ScenarioNames::SCENARIO6,
//         ScenarioNames::SCENARIO7,ScenarioNames::SCENARIO8,ScenarioNames::SCENARIO9,ScenarioNames::SCENARIO10,
//         ScenarioNames::SCENARIO11,ScenarioNames::SCENARIO12,ScenarioNames::SCENARIO13,ScenarioNames::SCENARIO14,ScenarioNames::SCENARIO15,
//         ScenarioNames::SCENARIO16,ScenarioNames::SCENARIO17,ScenarioNames::SCENARIO18,ScenarioNames::SCENARIO19,ScenarioNames::SCENARIO20
//     };
//     this->scenario_name = ScenarioNames::SCENARIO1; //default
//     this->probability = 0.05; // equal probability for each scenario
//     this->perturb = {
//         {ScenarioNames::SCENARIO1, 10.0}, 
//         {ScenarioNames::SCENARIO2, 20.0}, 
//         {ScenarioNames::SCENARIO3, 30.0},
//         {ScenarioNames::SCENARIO4, 40.0}, 
//         {ScenarioNames::SCENARIO5, 15}, 
//         {ScenarioNames::SCENARIO6, 25},
//         {ScenarioNames::SCENARIO7, 35}, 
//         {ScenarioNames::SCENARIO8, 45}, 
//         {ScenarioNames::SCENARIO9, 50},
//         {ScenarioNames::SCENARIO10, 5},
//         {ScenarioNames::SCENARIO11, 60}, 
//         {ScenarioNames::SCENARIO12, 65}, 
//         {ScenarioNames::SCENARIO13, 70},
//         {ScenarioNames::SCENARIO14, 75}, 
//         {ScenarioNames::SCENARIO15, 80}, 
//         {ScenarioNames::SCENARIO16, 90},
//         {ScenarioNames::SCENARIO17, 95}, 
//         {ScenarioNames::SCENARIO18, 55}, 
//         {ScenarioNames::SCENARIO19, 100},
//         {ScenarioNames::SCENARIO20, 105}
//         // {ScenarioNames::SCENARIO21, 0.05488135039273248},
//         // {ScenarioNames::SCENARIO22, 0.17151893663724196}, 
//         // {ScenarioNames::SCENARIO23, 0.2602763376071644},
//         // {ScenarioNames::SCENARIO24, 0.05488135039273248},
//         // {ScenarioNames::SCENARIO25, 0.17151893663724196}, 
//         // {ScenarioNames::SCENARIO26, 0.2602763376071644},
//         // {ScenarioNames::SCENARIO27, 0.05488135039273248},
//         // {ScenarioNames::SCENARIO28, 0.17151893663724196}, 
//         // {ScenarioNames::SCENARIO29, 0.2602763376071644},
//         // {ScenarioNames::SCENARIO30, 0.05488135039273248}


//     };
//     // this->first_stage_IX = {
//     //     mc::Interval(-1.1370129865018144, -1.1370129865018144),       
//     //     mc::Interval(1.5650443750201455, 1.5650443750201455),    
//     //     mc::Interval(-0.5986905053355759, -0.5986905053355759),
//     //     mc::Interval(1.3454091920547648, 1.3454091920547648),
//     //     mc::Interval(-0.23129676903369026, -0.23129676903369026),
//     //     mc::Interval(1.2019450763451505, 1.2019450763451505),
//     //     mc::Interval(0.052637388800998, 0.052637388800998),
//     //     mc::Interval(0.5559409367398721, 0.5559409367398721),
//     //     mc::Interval(0.6052476858593598, 0.6052476858593598),
//     //     mc::Interval(1.1, 1.1),
//     //     mc::Interval(0.8068686401687863, 0.8068686401687863),
//     //     mc::Interval(0.663161467338442, 0.663161467338442)
//     // };
//     // this->second_stage_IX = {
//     //     mc::Interval(5.137681424644748,5.137681424644748),       
//     //     mc::Interval(-4.95353950678979,-4.95353950678979),    
//     //     mc::Interval(2.6641201585529943,2.6641201585529943),      
//     //     mc::Interval(-2,-2), 
//     //     mc::Interval(1.885799773382196,1.885799773382196),

//     // };

//     this->first_stage_IX = {
//         mc::Interval(1e-5,16), //x5
//         mc::Interval(1e-5,16) //x6
//     };
        
//     this->second_stage_IX = {
//         mc::Interval(0, 1), //x1
//         mc::Interval(0, 1), //x2
//         mc::Interval(0, 1),  //x3
//         mc::Interval(0, 1) //x4

//     };
 
// };
// void Ex722Model::buildDAG() {
//     for (auto& scenario_name : this->scenario_names) {
//         int n_first_stage_vars = this->first_stage_IX.size();

//         // Loop over each scenario to build subproblem

//         const int nvars = n_first_stage_vars + this->second_stage_IX.size();

//         this->X[scenario_name].resize(nvars);


//         for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);

//         // scenario perturbation
//         double p = this->perturb[this->scenario_name];

//         mc::FFVar c1,c2,c3,c4,c5;
//         mc::FFVar nc1,nc2,nc3,nc4;

//         //.  0.09755988*x1*x5 + x1 == 1;
//         c1=0.09755988*this->X[scenario_name][2]*this->X[scenario_name][0]+this->X[scenario_name][2]-1;
//         nc1=-c1;

//         //.  0.0965842812*x2*x6 + x2 - x1 == 0;
//         c2=0.0965842812*this->X[scenario_name][3]*this->X[scenario_name][1]+this->X[scenario_name][3]-this->X[scenario_name][2];
//         nc2=-c2;

//         //.  0.0391908*x3*x5 + x3 + x1 == 1;
//         c3=0.0391908*this->X[scenario_name][4]*this->X[scenario_name][0]+this->X[scenario_name][4]+this->X[scenario_name][2]-1;
//         nc3=-c3;

//         // 0.03527172*x4*x6 + x4 - x1 + x2 - x3 == 0;
//         c4=0.03527172*this->X[scenario_name][5]*this->X[scenario_name][1]+this->X[scenario_name][5]-this->X[scenario_name][2]+this->X[scenario_name][3]-this->X[scenario_name][4];
//         nc4=-c4;

//         // x5**0.5 + x6**0.5 =L= 4;
//         c5=pow(this->X[scenario_name][0], 0.5)+pow(this->X[scenario_name][1], 0.5)-this->perturb[scenario_name];

//         mc::FFVar objective =-1000000*this->probability*this->X[scenario_name][5];
//         this->F[scenario_name]={objective,c1,c2,c3,c4,c5,nc1,nc2,nc3,nc4};
//     }
// }
// void Ex722Model::buildFullModelDAG(){
//     // for full model solve we will stay in scenario 1
//     int n_first_stage_vars = this->first_stage_IX.size();
//     int nvars = this->first_stage_IX.size()+this->second_stage_IX.size();
//     int n_second_stage_vars = this->second_stage_IX.size()/this->scenario_names.size(); // number of second stage variables per scenario
//     this->X[ScenarioNames::SCENARIO1].resize(this->first_stage_IX.size() +this->second_stage_IX.size());
    

//     for (int i = 0; i < n_first_stage_vars; ++i) this->X[ScenarioNames::SCENARIO1][i].set(&this->DAG[ScenarioNames::SCENARIO1]);
//     for (int s_idx=0; s_idx<this->scenario_names.size(); ++s_idx){
//         int second_stage_start_idx = n_first_stage_vars + s_idx * n_second_stage_vars;
//         for (int i = 0; i < n_second_stage_vars; ++i){
//             this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+i].set(&this->DAG[ScenarioNames::SCENARIO1]);
//         }
//     }
//     mc::FFVar objective=0;
//     for (int s_idx=0; s_idx<this->scenario_names.size(); ++s_idx){
//         int second_stage_start_idx = n_first_stage_vars + s_idx * n_second_stage_vars;

//         // scenario perturbation
//         double p = this->perturb[this->scenario_names[s_idx]];
//         mc::FFVar c1,c2,c3,c4,c5;
//         mc::FFVar nc1,nc2,nc3,nc4;

//         //.  0.09755988*x1*x5 + x1 == 1;
//         c1=0.09755988*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]*this->X[ScenarioNames::SCENARIO1][0]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]-1;
//         nc1=-c1;

//         //.  0.0965842812*x2*x6 + x2 - x1 == 0;
//         c2=0.0965842812*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]*this->X[ScenarioNames::SCENARIO1][1]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]-this->X[ScenarioNames::SCENARIO1][second_stage_start_idx];
//         nc2=-c2;

//         //.  0.0391908*x3*x5 + x3 + x1 == 1;
//         c3=0.0391908*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+2]*this->X[ScenarioNames::SCENARIO1][0]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+2]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]-1;
//         nc3=-c3;

//         // 0.03527172*x4*x6 + x4 - x1 + x2 - x3 == 0;
//         c4=0.03527172*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+3]*this->X[ScenarioNames::SCENARIO1][1]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+3]-this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]-this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+2];
//         nc4=-c4;

//         // x5**0.5 + x6**0.5 =L= 4;
//         c5=pow(this->X[ScenarioNames::SCENARIO1][0], 0.5)+pow(this->X[ScenarioNames::SCENARIO1][1], 0.5)-this->perturb[this->scenario_names[s_idx]];

//         objective +=-1000000*this->probability*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+3];

//         std::vector<mc::FFVar> scenario_constraints = {c1,c2,c3,c4,nc1,nc2,nc3,nc4,c5};
//         this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].end(), scenario_constraints.begin(), scenario_constraints.end());
//     }
//     this->F[ScenarioNames::SCENARIO1].insert(this->F[ScenarioNames::SCENARIO1].begin(), objective);
//     this->full_model_built = true;
// }
// Ipopt::SmartPtr<STModel> Ex722Model::clone(){
//     Ipopt::SmartPtr<Ex722Model> p = new Ex722Model();

//     p->scenario_name=this->scenario_name;
//     p->first_stage_IX=this->first_stage_IX;
//     p->second_stage_IX=this->second_stage_IX;
//     p->perturb=this->perturb;
//     p->perturb_coeffs=this->perturb_coeffs;
//     p->scenario_names=this->scenario_names;
//     p->probability=this->probability;
//     p->first_stage_map=this->first_stage_map;
//     p->second_stage_map=this->second_stage_map;
//     p->clearDAG(); // clear the DAG of the cloned model
//     if (this->full_model_built) {
//         p->buildFullModelDAG();
//     } else {
//         p->buildDAG();
//     }
//     return p;
// }
#include "Ex722.h"

Ex722Model::Ex722Model(BranchingStrategy branching_strategy):STModel() {

    this->branching_strategy = branching_strategy;
    this->scenario_names = { ScenarioNames::SCENARIO1, ScenarioNames::SCENARIO2,ScenarioNames::SCENARIO3,
        ScenarioNames::SCENARIO4, ScenarioNames::SCENARIO5
    };
    this->scenario_name = ScenarioNames::SCENARIO1; //default
    this->probability = 0.2; // equal probability for each scenario

    // feed-type driver -> e5's bound directly, same values as before
    this->perturb = {
        {ScenarioNames::SCENARIO1, 10.0},
        {ScenarioNames::SCENARIO2, 20.0},
        {ScenarioNames::SCENARIO3, 30.0},
        {ScenarioNames::SCENARIO4, 40.0},
        {ScenarioNames::SCENARIO5, 15}
    };

    // ------------------------------------------------------------------
    // Scenario data matching the Python const_model(): explicit lists,
    // one value per scenario, for the underlying uncertain quantities.
    // temp_factor scales the bilinear coupling coefficients (a1-a4);
    // price_shock scales the second-stage recourse cost. Every
    // scenario-dependent coefficient below is derived from these lists,
    // in the same order as this->scenario_names.
    // ------------------------------------------------------------------
    std::vector<double> temp_factor = {
        0.85, 0.90, 0.95, 1.00, 1.05
    };
    std::vector<double> price_shock = {
        1.00, 1.05, 0.95, 1.10, 0.90
    };

    for (size_t i = 0; i < this->scenario_names.size(); ++i) {
        ScenarioNames sn = this->scenario_names[i];
        double tf = temp_factor[i];

        // {a1, a2, a3, a4} for this scenario, matching Python's
        // a1_vals/a2_vals/a3_vals/a4_vals
        this->perturb_coeffs[sn] = {
            0.09755988 * tf,             // a1
            0.0965842812 * std::sqrt(tf),// a2
            0.0391908 * tf,              // a3
            0.03527172 * std::sqrt(tf)   // a4
        };

        this->recourse_cost[sn] = 1000.0 * price_shock[i];
    }

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
    for (auto& scenario_name : this->scenario_names) {
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
    int n_second_stage_vars = this->second_stage_IX.size(); // per-scenario count (2: x4, x5)
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