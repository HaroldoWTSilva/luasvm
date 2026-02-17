local svm = require "svm"

math.randomseed(2345678)
local test = require "test"
local kernels = require "kernel"

local C = 1

local fitfunc = function(tset)
	return svm.fit(tset, kernels.get_poly(0.5, 4), C )
end

test.run(test.get_generator_quadcircle, 8192, fitfunc)
