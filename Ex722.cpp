#include "Ex722.h"

Ex722Model::Ex722Model(BranchingStrategy branching_strategy):STModel() {

    this->branching_strategy = branching_strategy;
    this->scenario_names = { ScenarioNames::SCENARIO1, ScenarioNames::SCENARIO2,ScenarioNames::SCENARIO3,
     ScenarioNames::SCENARIO4, ScenarioNames::SCENARIO5,ScenarioNames::SCENARIO6,
     ScenarioNames::SCENARIO7,ScenarioNames::SCENARIO8,ScenarioNames::SCENARIO9,ScenarioNames::SCENARIO10
    };
    this->scenario_name = ScenarioNames::SCENARIO1; //default
    this->probability = 0.1; // equal probability for each scenario
    this->perturb = {
        {ScenarioNames::SCENARIO1, 4}, 
        {ScenarioNames::SCENARIO2, 3}, 
        {ScenarioNames::SCENARIO3, 4},
        {ScenarioNames::SCENARIO4, 3}, 
        {ScenarioNames::SCENARIO5, 4}, 
        {ScenarioNames::SCENARIO6, 3},
        {ScenarioNames::SCENARIO7, 4}, 
        {ScenarioNames::SCENARIO8, 3}, 
        {ScenarioNames::SCENARIO9, 4},
        {ScenarioNames::SCENARIO10, 3}
        // {ScenarioNames::SCENARIO11, -1}, 
        // {ScenarioNames::SCENARIO12, -1}, 
        // {ScenarioNames::SCENARIO13, 0.2602763376071644},
        // {ScenarioNames::SCENARIO14, 0.05488135039273248}, 
        // {ScenarioNames::SCENARIO15, 0.17151893663724196}, 
        // {ScenarioNames::SCENARIO16, 0.2602763376071644},
        // {ScenarioNames::SCENARIO17, 0.05488135039273248}, 
        // {ScenarioNames::SCENARIO18, 0.17151893663724196}, 
        // {ScenarioNames::SCENARIO19, 0.2602763376071644},
        // {ScenarioNames::SCENARIO20, 0.05488135039273248}
        // {ScenarioNames::SCENARIO21, 0.05488135039273248},
        // {ScenarioNames::SCENARIO22, 0.17151893663724196}, 
        // {ScenarioNames::SCENARIO23, 0.2602763376071644},
        // {ScenarioNames::SCENARIO24, 0.05488135039273248},
        // {ScenarioNames::SCENARIO25, 0.17151893663724196}, 
        // {ScenarioNames::SCENARIO26, 0.2602763376071644},
        // {ScenarioNames::SCENARIO27, 0.05488135039273248},
        // {ScenarioNames::SCENARIO28, 0.17151893663724196}, 
        // {ScenarioNames::SCENARIO29, 0.2602763376071644},
        // {ScenarioNames::SCENARIO30, 0.05488135039273248}


    };
    // this->first_stage_IX = {
    //     mc::Interval(-1.1370129865018144, -1.1370129865018144),       
    //     mc::Interval(1.5650443750201455, 1.5650443750201455),    
    //     mc::Interval(-0.5986905053355759, -0.5986905053355759),
    //     mc::Interval(1.3454091920547648, 1.3454091920547648),
    //     mc::Interval(-0.23129676903369026, -0.23129676903369026),
    //     mc::Interval(1.2019450763451505, 1.2019450763451505),
    //     mc::Interval(0.052637388800998, 0.052637388800998),
    //     mc::Interval(0.5559409367398721, 0.5559409367398721),
    //     mc::Interval(0.6052476858593598, 0.6052476858593598),
    //     mc::Interval(1.1, 1.1),
    //     mc::Interval(0.8068686401687863, 0.8068686401687863),
    //     mc::Interval(0.663161467338442, 0.663161467338442)
    // };
    // this->second_stage_IX = {
    //     mc::Interval(5.137681424644748,5.137681424644748),       
    //     mc::Interval(-4.95353950678979,-4.95353950678979),    
    //     mc::Interval(2.6641201585529943,2.6641201585529943),      
    //     mc::Interval(-2,-2), 
    //     mc::Interval(1.885799773382196,1.885799773382196),

    // };

    this->first_stage_IX = {
        mc::Interval(0, 1),
        mc::Interval(0, 1),
        mc::Interval(0, 1),
        mc::Interval(0, 1)
    };
        
    this->second_stage_IX = {
        mc::Interval(1e-5,16),
        mc::Interval(1e-5,16)

    };
 
};
void Ex722Model::buildDAG() {
    for (auto& scenario_name : this->scenario_names) {
        int n_first_stage_vars = this->first_stage_IX.size();

        // Loop over each scenario to build subproblem

        const int nvars = n_first_stage_vars + this->second_stage_IX.size();

        this->X[scenario_name].resize(nvars);


        for (int i = 0; i < nvars; ++i) this->X[scenario_name][i].set(&this->DAG[scenario_name]);

        // scenario perturbation
        double p = this->perturb[this->scenario_name];

        mc::FFVar c1,c2,c3,c4,c5;
        mc::FFVar nc1,nc2,nc3,nc4;

        //.  0.09755988*x1*x5 + x1 == 1;
        c1=0.09755988*this->X[scenario_name][0]*this->X[scenario_name][4]+this->X[scenario_name][0]-1;
        nc1=-c1;

        //.  0.0965842812*x2*x6 + x2 - x1 == 0;
        c2=0.0965842812*this->X[scenario_name][1]*this->X[scenario_name][5]+this->X[scenario_name][1]-this->X[scenario_name][0];
        nc2=-c2;

        //.  0.0391908*x3*x5 + x3 + x1 == 1;
        c3=0.0391908*this->X[scenario_name][2]*this->X[scenario_name][4]+this->X[scenario_name][2]+this->X[scenario_name][0]-1;
        nc3=-c3;

        // 0.03527172*x4*x6 + x4 - x1 + x2 - x3 == 0;
        c4=0.03527172*this->X[scenario_name][3]*this->X[scenario_name][5]+this->X[scenario_name][3]-this->X[scenario_name][0]+this->X[scenario_name][1]-this->X[scenario_name][2];
        nc4=-c4;

        // x5**0.5 + x6**0.5 =L= 4;
        c5=pow(this->X[scenario_name][4], 0.5)+pow(this->X[scenario_name][5], 0.5)-this->perturb[scenario_name];

        mc::FFVar objective =-10000*this->probability*this->X[scenario_name][3];
        this->F[scenario_name]={objective,c1,c2,c3,c4,c5,nc1,nc2,nc3,nc4};
    }
}
void Ex722Model::buildFullModelDAG(){
    // for full model solve we will stay in scenario 1
    int n_first_stage_vars = this->first_stage_IX.size();
    int nvars = this->first_stage_IX.size()+this->second_stage_IX.size();
    int n_second_stage_vars = this->second_stage_IX.size()/this->scenario_names.size(); // number of second stage variables per scenario
    this->X[ScenarioNames::SCENARIO1].resize(this->first_stage_IX.size() +this->second_stage_IX.size());
    

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

        // scenario perturbation
        double p = this->perturb[this->scenario_names[s_idx]];
        mc::FFVar c1,c2,c3,c4,c5;
        mc::FFVar nc1,nc2,nc3,nc4;

        //.  0.09755988*x1*x5 + x1 == 1;
        c1=0.09755988*this->X[ScenarioNames::SCENARIO1][0]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][0]-1;
        nc1=-c1;

        //.  0.0965842812*x2*x6 + x2 - x1 == 0;
        c2=0.0965842812*this->X[ScenarioNames::SCENARIO1][1]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]+this->X[ScenarioNames::SCENARIO1][1]-this->X[ScenarioNames::SCENARIO1][0];
        nc2=-c2;

        //.  0.0391908*x3*x5 + x3 + x1 == 1;
        c3=0.0391908*this->X[ScenarioNames::SCENARIO1][2]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx]+this->X[ScenarioNames::SCENARIO1][2]+this->X[ScenarioNames::SCENARIO1][0]-1;
        nc3=-c3;

        // 0.03527172*x4*x6 + x4 - x1 + x2 - x3 == 0;
        c4=0.03527172*this->X[ScenarioNames::SCENARIO1][3]*this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1]+this->X[ScenarioNames::SCENARIO1][3]-this->X[ScenarioNames::SCENARIO1][0]+this->X[ScenarioNames::SCENARIO1][1]-this->X[ScenarioNames::SCENARIO1][2];
        nc4=-c4;

        // x5**0.5 + x6**0.5 =L= 4;
        c5=pow(this->X[ScenarioNames::SCENARIO1][second_stage_start_idx], 0.5)+pow(this->X[ScenarioNames::SCENARIO1][second_stage_start_idx+1], 0.5)-this->perturb[this->scenario_names[s_idx]];

        objective +=-10000*this->probability*this->X[ScenarioNames::SCENARIO1][3];

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
