function Population = MygacreationuniformForSTWE(GenomeLength,FitnessFcn,options)
% add stwe_setup into the file gacreationuniform, 20100131 by TYB

%GACREATIONUNIFORM Creates the initial population for genetic algorithm.
%   POP = GACREATIONUNIFORM(NVARS,FITNESSFCN,OPTIONS) Creates the
%   initial population that GA will then evolve into a solution.
%
%   Population size can be a vector of separate populations.
%   Here, we are only interested in the total number.
%
%   Example:
%     options = gaoptimset('PopulationType','bitString');
%            NVARS = 1; FitnessFcn = @ackleyfcn;
%
%     pop = gacreationuniform(NVARS,FitnessFcn,options);
%
%   pop will be a 20-by-1 logical column vector.  Note that the
%   default Population Size in GAOPTIMSET is 20.

%   Copyright 2003-2006 The MathWorks, Inc.
%   $Revision: 1.7.4.4 $  $Date: 2006/06/20 20:06:41 $
%

if strcmpi(options.PopulationType,'custom')
    error('gads:GACREATIONUNIFORM:unknownPopulationType','Can not create initial population for ''%s'' population type.',options.PopulationType);
end

stwe_setup();  %% add by TYB 20100131

totalPopulation = sum(options.PopulationSize);
initPopProvided = size(options.InitialPopulation,1);
individualsToCreate = totalPopulation - initPopProvided;

if strcmpi(options.PopulationType,'doubleVector')
    % Iitialize Population to be created
    Population = zeros(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    % problemtype is either 'unconstrained', 'boundconstraints', or
    % 'linearconstraints'. Nonlinear constrained algorithm 'ALGA' does not
    % create or use initial population of its own. It calls sub-problem
    % solvers (galincon/gaunc)
    if isfield(options,'LinearConstr')
        problemtype = options.LinearConstr.type;
    else
        problemtype = 'unconstrained';
    end
    % This function knows how to create initial population for
    % unconstrained and boundconstrained cases but does not know in
    % linearconstrained case.
    if ~strcmp(problemtype,'linearconstraints')
        range = options.PopInitRange;
        lowerBound = range(1,:);
        span = range(2,:) - lowerBound;
        Population(initPopProvided+1:end,:) = repmat(lowerBound,individualsToCreate,1) + ...
            repmat(span,individualsToCreate,1) .* rand(individualsToCreate,GenomeLength);
    else
        Population(initPopProvided+1:end,:) = []; 
    end
elseif strcmpi(options.PopulationType,'bitString')
    % Iitialize Population to be created
    Population = true(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
    Population(initPopProvided+1:end,:) = logical(0.5 > rand(individualsToCreate,GenomeLength));
end

if all(isnan(Population))
    error('gads:GACREATIONUNIFORM:populationIsNaN',['Initial population contains NaN;','OPTIONS.PopInitRange is possibly too big.']);
elseif all(isinf(Population))
    error('gads:GACREATIONUNIFORM:populationIsInf',['Initial population contains Inf;','OPTIONS.PopInitRange is possibly too big.']);
end

