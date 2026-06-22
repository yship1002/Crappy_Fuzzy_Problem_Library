#include "STModel.h"

class EDUnits_nocp : public STModel {
    public:
        EDUnits_nocp(BranchingStrategy branchingStrategy);
        EDUnits_nocp(const EDUnits_nocp& other) = default;
        EDUnits_nocp() = default;

        double probability;
        Ipopt::SmartPtr<STModel> clone() override;
        void buildDAG() override;
        void buildFullModelDAG() override;

        template<class Archive>
        void serialize(Archive& ar) {
            CEREAL_NVP(first_stage_IX), CEREAL_NVP(second_stage_IX);
        }
};