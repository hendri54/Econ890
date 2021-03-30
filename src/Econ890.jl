module Econ890

using DocStringExtensions, Interpolations, Roots

export AbstractUtility, UtilityCRRA, UtilityLog
export euler_dev, utility, marg_utility, inv_utility, inv_marg_utility, c_growth

# Policy function iteration
export Model, init_test_model, euler_dev, solve

include("utility.jl");
include("vfi.jl");

end # module
