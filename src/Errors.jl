# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

using PowerDynBase: PowerDynamicsError

struct GridSolutionError <: PowerDynamicsError
    msg::String
end
