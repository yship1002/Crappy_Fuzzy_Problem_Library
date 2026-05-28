#include "STModel.h"

class EDUnits : public STModel {
    public:
        EDUnits(BranchingStrategy branchingStrategy);
        EDUnits(const EDUnits& other) = default;
        EDUnits() = default;

        double probability;
        Ipopt::SmartPtr<STModel> clone() override;
        void buildDAG() override;
        void buildFullModelDAG() override;

        template<class Archive>
        void serialize(Archive& ar) {
            CEREAL_NVP(first_stage_IX), CEREAL_NVP(second_stage_IX);
        }
};