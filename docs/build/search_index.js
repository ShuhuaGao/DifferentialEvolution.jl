var documenterSearchIndex = {"docs":
[{"location":"common/#Common-interfaces-of-DE","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"","category":"section"},{"location":"common/#Core-DE","page":"Common interfaces of DE","title":"Core DE","text":"","category":"section"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"The core data structure of differential evolution is a struct that inherits the abstract type Evolver or, more typically, the DEEvolver, which itself inherits Evolver and  represents a single-population DE. ","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"Evolver\r\nOptions\r\nDEEvolver\r\nDEOptions","category":"page"},{"location":"common/#DifferentialEvolution.Evolver","page":"Common interfaces of DE","title":"DifferentialEvolution.Evolver","text":"Evolver\n\nAbstract type of all DE core.\n\n\n\n\n\n","category":"type"},{"location":"common/#DifferentialEvolution.Options","page":"Common interfaces of DE","title":"DifferentialEvolution.Options","text":"Options\n\nAbstract type of all DE options.\n\n\n\n\n\n","category":"type"},{"location":"common/#DifferentialEvolution.DEEvolver","page":"Common interfaces of DE","title":"DifferentialEvolution.DEEvolver","text":"DEEvolver <: Evolver\n\nAn abstract normal (single-population) DE evolution type \n\n\n\n\n\n","category":"type"},{"location":"common/#DifferentialEvolution.DEOptions","page":"Common interfaces of DE","title":"DifferentialEvolution.DEOptions","text":"DEOptions <: Options\n\nAbstract type of options for a DE, particually, a single-population one. See also DEEvolver.\n\n\n\n\n\n","category":"type"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"Each concrete type of DEEvolver should support the five default interfaces listed below. If the five corresponding fields are defined, then nothing extra needs to be done.","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"pop(e::DEEvolver) = e.pop\r\nfitness(e::DEEvolver) = e.fitness\r\nsenses(e::DEEvolver) = e.senses\r\nbounds(e::DEEvolver) = e.bounds\r\noptions(e::DEEvolver) = e.options","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"Thus, the five fundamental fields listed above can be retrieved directly or using a method with the same name of the field.","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"pop\r\nfitness\r\nsenses\r\nbounds\r\noptions","category":"page"},{"location":"common/#DifferentialEvolution.pop","page":"Common interfaces of DE","title":"DifferentialEvolution.pop","text":"pop(e::DEEvolver)\n\nGet the population of a DE, which is typically a matrix.\n\n\n\n\n\n","category":"function"},{"location":"common/#DifferentialEvolution.fitness","page":"Common interfaces of DE","title":"DifferentialEvolution.fitness","text":"fitness(e::DEEvolver)\n\nGet the fitness of the population, which is generally a vector for single-objective DE or a matrix for multi-objective DE.\n\n\n\n\n\n","category":"function"},{"location":"common/#DifferentialEvolution.senses","page":"Common interfaces of DE","title":"DifferentialEvolution.senses","text":"senses(e::DEEvolver)\n\nGet the senses of optimization. In traditional evolutionary algorithms, the fitness of an individual is going to be maximized. If minimization is required instead, then set a negative sense.\n\n\n\n\n\n","category":"function"},{"location":"common/#DifferentialEvolution.bounds","page":"Common interfaces of DE","title":"DifferentialEvolution.bounds","text":"bounds(e::DEEvolver)\n\nGet the bound constraints,typically a n x 2 matrix for n variables. \n\n\n\n\n\n","category":"function"},{"location":"common/#DifferentialEvolution.options","page":"Common interfaces of DE","title":"DifferentialEvolution.options","text":"options(e::DEEvolver)\n\nGet the options of this DE.\n\n\n\n\n\n","category":"function"},{"location":"common/#DE-options","page":"Common interfaces of DE","title":"DE options","text":"","category":"section"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"Options, i.e., control parameters, of DE are encapsulated in a struct as a subtype of DEOptions (for the most common single-population DE) or a subtype of Options. The options of a given DE is accessed by the method options.","category":"page"},{"location":"common/#Other-utility-interfaces","page":"Common interfaces of DE","title":"Other utility interfaces","text":"","category":"section"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"The following three methods have default implementation for DEEvolver.","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"nindividuals(e::DEEvolver) \r\nndims(de::DEEvolver)\r\nnobjectives(e::DEEvolver)","category":"page"},{"location":"common/#DifferentialEvolution.nindividuals-Tuple{DEEvolver}","page":"Common interfaces of DE","title":"DifferentialEvolution.nindividuals","text":"nindividuals(e::DEEvolver)\n\nNumber of individuals in the population.\n\n\n\n\n\n","category":"method"},{"location":"common/#Base.ndims-Tuple{DEEvolver}","page":"Common interfaces of DE","title":"Base.ndims","text":"ndims(de::DEEvolver)\n\nNumber of variables (dimension) of the DE.\n\n\n\n\n\n","category":"method"},{"location":"common/#DifferentialEvolution.nobjectives-Tuple{DEEvolver}","page":"Common interfaces of DE","title":"DifferentialEvolution.nobjectives","text":"nobjectives(e::Evolver)\n\nNumber of objectives.\n\n\n\n\n\n","category":"method"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"No default implementation is provided for the two methods below, since a correct implementation depends on details of the concrete DE type.","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"best_individual(e::Evolver; all=false)\r\nbest_fitness(e::Evolver)","category":"page"},{"location":"common/#DifferentialEvolution.best_individual-Tuple{Evolver}","page":"Common interfaces of DE","title":"DifferentialEvolution.best_individual","text":"best_individual(e::Evolver; all=false)\n\nGet individuals that have the best fitness. If all is true, then all such individuals are returned; otherwise, only one of them is fetched.\n\n\n\n\n\n","category":"method"},{"location":"common/#DifferentialEvolution.best_fitness-Tuple{Evolver}","page":"Common interfaces of DE","title":"DifferentialEvolution.best_fitness","text":"best_fitness(e::Evolver)\n\nGet the best fitness among the individuals of the evolver.\n\n\n\n\n\n","category":"method"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"It is highly recommended that a concrete DE type implements the two best methods as outlined above.","category":"page"},{"location":"common/","page":"Common interfaces of DE","title":"Common interfaces of DE","text":"Currently, DifferentialEvolution.jl provides a standard implementation with the types StdDE and StdOptions.","category":"page"},{"location":"stdDE/#Standard-DE","page":"Standard DE","title":"Standard DE","text":"","category":"section"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"The standard single-objective DE types and evolution algortithms. ","category":"page"},{"location":"stdDE/#Types","page":"Standard DE","title":"Types","text":"","category":"section"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"StdOptions\r\nStdDE\r\nStdSoDE\r\nStdMoDE","category":"page"},{"location":"stdDE/#DifferentialEvolution.StdOptions","page":"Standard DE","title":"DifferentialEvolution.StdOptions","text":"StdOptions{T<:AF}\n\nStandard options with two parameters, F (scaling) and Cr (crossover rate). This type can be constructed with either positional arguments or keyword arguments.\n\n\n\n\n\n","category":"type"},{"location":"stdDE/#DifferentialEvolution.StdDE","page":"Standard DE","title":"DifferentialEvolution.StdDE","text":"StdDE{TP<:AM{<:AF}, TF<:AVoM{<:Real}, TS<:ASoV{<:Integer}, TO<:DEOptions} <: DEEvolver\n\nThe standard DE type with the following fields:\n\npop::TP: a matrix-type population, each column of which is an individual (vector) in DE\nfitness::TF: fitness of each individual. For single-objective DE, fitness is a vector-like    structure; for multi-objective DE, fitness is of a matrix type, whose column contains values    of multiple objectives for each individual.\nsenses::TS: integer for a single objective or vector of integers for multiple objectives.    An individual's strength is defined by fitness times senses. Only values 1 (maximization)    and -1 (minimization) are allowed. To perform minimization of an    objective, set a negative sense, i.e., -1.\nbounds::TP: a d-by-2 matrix for a d-dim DE, each row containing the lower and upper boundaries   of a variable.\noptions::TO: options which usually contain the control parameters. See also StdOptions.\n\nIn practice, more convenient types are usually used like StdSoDE for single-objective DE  and StdMoDE for multiobjective DE.\n\n\n\n\n\n","category":"type"},{"location":"stdDE/#DifferentialEvolution.StdSoDE","page":"Standard DE","title":"DifferentialEvolution.StdSoDE","text":"StdSoDE{TP<:AM{<:AF}, TF<:AV{<:Real}, TW<:Real, TO<:DEOptions}\n\nStandard single-objective DE.\n\n\n\n\n\n","category":"type"},{"location":"stdDE/#DifferentialEvolution.StdMoDE","page":"Standard DE","title":"DifferentialEvolution.StdMoDE","text":"StdMoDE{TP<:AM{<:AF}, TF<:AM{<:Real}, TW<:AV{<:Real}, TO<:DEOptions}\n\nStandard multi-objective DE.\n\n\n\n\n\n","category":"type"},{"location":"stdDE/#Constructor","page":"Standard DE","title":"Constructor","text":"","category":"section"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"The three standard DE types StdDE, StdSoDE, and StdMoDE can all be built with their inner constructors, but a more convenient way is to employ the dedicated outer constructor. Please check Julia documentation on Constructors for more details.","category":"page"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"CurrentModule = DifferentialEvolution","category":"page"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; \r\n        senses::ASoV{<:Integer}) where {TF<:Real}\r\n\r\nStdDE(nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; senses::ASoV{<:Integer})","category":"page"},{"location":"stdDE/#DifferentialEvolution.StdDE-Union{Tuple{TF}, Tuple{Type{TF},Int64,AbstractArray{var\"#s17\",2} where var\"#s17\"<:AbstractFloat,DEOptions}} where TF<:Real","page":"Standard DE","title":"DifferentialEvolution.StdDE","text":"StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AF}, opt::DEOptions; weights::ASoV{<:Real}) where {TF<:Real}\n\nCreate a standard DE evolver with a population of nindividuals individuals and options specified by  opt. The length of each individual (also known as a vector) is equal to the number of rows of bounds. That is, each row of bounds denotes the range l u of a variable, where l and u denote  the lower and upper bound, respectively. \n\nTF is the data type of fitness. The type of individuals (solution vectors) is the same as the type of bounds.\n\nThe size of senses indicates the number of objectives. If the supplied senses is a scalar, then a single-objective DE StdSoDE instance is built; otherwise, a multi-objective StdMoDE is constructed. \n\n\n\n\n\n","category":"method"},{"location":"stdDE/#DifferentialEvolution.StdDE-Tuple{Int64,AbstractArray{var\"#s17\",2} where var\"#s17\"<:AbstractFloat,DEOptions}","page":"Standard DE","title":"DifferentialEvolution.StdDE","text":"StdDE(nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions; senses::ASoV{<:Integer})\n\nConstruct a standard DE whose vectors are of default type Float64. See StdDE(::Type{TF}, nindividuals::Int, bounds::AM{<:AbstractFloat}, opt::DEOptions;          senses::ASoV{<:Integer}) where {TF<:Real} for details.\n\n\n\n\n\n","category":"method"},{"location":"stdDE/#Evolution","page":"Standard DE","title":"Evolution","text":"","category":"section"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"The basic single-objective DE evolution algorithm is implemented.","category":"page"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"evolve!","category":"page"},{"location":"stdDE/#DifferentialEvolution.evolve!","page":"Standard DE","title":"DifferentialEvolution.evolve!","text":"evolve!(de::StdDES, evaluator::Function, ngen; mut=mutate_rand_1, cx=crossover_binomial,\n    bc=:bounce_back)\n\nStandard single-objective differential evolution. \n\nEach individual's fitness is evaluated by the given evaluator and the evolution lasts for ngen iterations. Keyword arguments:\n\nmut=mutate_rand_1!: mutation method\ncx=crossover_binomial!: crossover method, currently support\nbc=:bounce_back: boundary constraints handling method， currently support\n:bounce_back\nstats=[:best, :avg]: fitness statistics to display in each generation.    Valid entries are :best, :avg(average), :max, :min, :std(standard deviation).\nstats_fmt::FormatSpec=FormatSpec(\"10.5e\"): format specification of each statistic term for display.   See Formatting.jl for details.\n\n\n\n\n\n","category":"function"},{"location":"stdDE/#Ranking","page":"Standard DE","title":"Ranking","text":"","category":"section"},{"location":"stdDE/","page":"Standard DE","title":"Standard DE","text":"best_individual(de::StdSoDE; all=false)\r\nbest_fitness(de::StdSoDE)\r\nBase.sort!(de::StdSoDE; rev::Bool=false)","category":"page"},{"location":"stdDE/#DifferentialEvolution.best_individual-Tuple{StdDE{TP,TF,TW,TO} where TO<:DEOptions where TW<:Real where TF<:(AbstractArray{var\"#s32\",1} where var\"#s32\"<:Real) where TP<:(AbstractArray{var\"#s33\",2} where var\"#s33\"<:AbstractFloat)}","page":"Standard DE","title":"DifferentialEvolution.best_individual","text":"best_individual(de::StdSoDE; all=false)\n\nGet the best individual (or individuals) from the single-objective DE population. See also best_individual(e::Evolver; all=false) for more details.\n\n\n\n\n\n","category":"method"},{"location":"stdDE/#DifferentialEvolution.best_fitness-Tuple{StdDE{TP,TF,TW,TO} where TO<:DEOptions where TW<:Real where TF<:(AbstractArray{var\"#s32\",1} where var\"#s32\"<:Real) where TP<:(AbstractArray{var\"#s33\",2} where var\"#s33\"<:AbstractFloat)}","page":"Standard DE","title":"DifferentialEvolution.best_fitness","text":"best_fitness(de::StdSoDE)\n\nGet the best fitness from the single-objective DE population.\n\n\n\n\n\n","category":"method"},{"location":"stdDE/#Base.sort!-Tuple{StdDE{TP,TF,TW,TO} where TO<:DEOptions where TW<:Real where TF<:(AbstractArray{var\"#s32\",1} where var\"#s32\"<:Real) where TP<:(AbstractArray{var\"#s33\",2} where var\"#s33\"<:AbstractFloat)}","page":"Standard DE","title":"Base.sort!","text":"sort!(de::StdDES; rev::Bool=false)\n\nSort a single-objective DE population in place according to the fitness. The best individual appears at first if rev is false (the default value).\n\n\n\n\n\n","category":"method"},{"location":"mutation/#Mutation","page":"Mutation","title":"Mutation","text":"","category":"section"},{"location":"mutation/","page":"Mutation","title":"Mutation","text":"Common mutation operators in DE are provided. We recommend in-place operators for efficiency purpose. That is, a dummy donor vector is the first argument, which is then filled inside the method.","category":"page"},{"location":"mutation/#DE/rand/1","page":"Mutation","title":"DE/rand/1","text":"","category":"section"},{"location":"mutation/","page":"Mutation","title":"Mutation","text":"mutate_rand_1!","category":"page"},{"location":"mutation/#DifferentialEvolution.mutate_rand_1!","page":"Mutation","title":"DifferentialEvolution.mutate_rand_1!","text":"mutate_rand_1!(V::AV{T}, i::Integer, de::StdDE{<:AM{T}}) where T<:AF\n\nDE/rand/1 mutation in place, which produces a donor vector V to crossover with the i-th target vector (parent).\n\n\n\n\n\nmutate_rand_1!(V::AV{T}, i::Integer, pop::AM{T}, F::AF, indices::AV{<:Integer}) where T<:AF\n\nA lower-level implementation. V is the donor vector to be filled, i is the index of the target  vector in the population pop whose size is n, F is the scale factor, and indices is an  array of length n that contains intergers from 1 to n (in an arbitrary order).\n\n\n\n\n\n","category":"function"},{"location":"crossover/#Crossover","page":"Crossover","title":"Crossover","text":"","category":"section"},{"location":"crossover/","page":"Crossover","title":"Crossover","text":"Crossover (also known as mating or recombination) operations in DE. Note that we recommend  in-place crossover operators: supposing two vectors X and V are recombined into another  vector U, the vector U should be provided as the first argument and be filled inside the method.","category":"page"},{"location":"crossover/#Binomial-crossover","page":"Crossover","title":"Binomial crossover","text":"","category":"section"},{"location":"crossover/","page":"Crossover","title":"Crossover","text":"crossover_binomial!","category":"page"},{"location":"crossover/#DifferentialEvolution.crossover_binomial!","page":"Crossover","title":"DifferentialEvolution.crossover_binomial!","text":"crossover_binomial!(U::AV{T}, i::Integer, de::StdDE{<:AM{T}}, V::AV{T}) where T<:AF\n\nRecombine the i-th vector of de and the donor vector V into a trial vector U. The crossover rate is specied by de.options. \n\nSee also StdDE.\n\n\n\n\n\ncrossover_binomial!(U::AV{T}, X::AV{T}, V::AV{T}, Cr::AF) where T<:AF\n\nRecombine two vectors X and V into U subject to the rate Cr.\n\n\n\n\n\n","category":"function"},{"location":"#Documentation-of-DifferentialEvolution.jl","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"","category":"section"},{"location":"#Type-aliases","page":"Documentation of DifferentialEvolution.jl","title":"Type aliases","text":"","category":"section"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"To avoid writing lengthy type names, DifferentialEvolution.jl defines the following abbreviations.","category":"page"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"const AbstractScalarOrVec{T} = Union{T, AbstractVector{T}} where T<:Real\r\nconst AV = AbstractVector\r\nconst AM = AbstractMatrix\r\nconst AVoM = AbstractVecOrMat\r\nconst ASoV = AbstractScalarOrVec\r\nconst AF = AbstractFloat","category":"page"},{"location":"#Contents","page":"Documentation of DifferentialEvolution.jl","title":"Contents","text":"","category":"section"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"Common interfaces of DE","category":"page"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"Mutation","category":"page"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"Crossover","category":"page"},{"location":"","page":"Documentation of DifferentialEvolution.jl","title":"Documentation of DifferentialEvolution.jl","text":"Standard DE","category":"page"}]
}
