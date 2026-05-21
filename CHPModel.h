#include "STModel.h"

class CHPModel : public STModel {
public:
    CHPModel(BranchingStrategy branching_strategy);
    CHPModel(const CHPModel& other) = default;
    CHPModel() = default;

    double probability;
    double invest_coeff; // 1e-6 * (af + c_m) * c_ref

    Ipopt::SmartPtr<STModel> clone() override;
    void buildDAG() override;
    void buildFullModelDAG() override;

    template<class Archive>
    void serialize(Archive& ar) {
        CEREAL_NVP(first_stage_IX), CEREAL_NVP(second_stage_IX);
    }
};
