using CSV,DelimitedFiles,DataFrames,Random,Distributions

include("species_level_network.jl")
include("./SLN_maker.jl")
include("./r_no_prey.jl")

P = CSV.read("Anisian_guild_parameters_11_12_19.csv",DataFrame)
A = readdlm("Tethys_B_matrix_11_08_19.csv")

outfile1 = ARGS[1]
f1 = open(outfile1, "w")

γ = 2.5

begin
    # calculate metanetwork diversity
    #no. of guilds
    no_guilds = size(P,1)
    println("\nNo. of guilds, G = ", no_guilds)
    #calculate number of species
    S = sum(P[:,:G])
    no_species = S[1]
    println("No. of species, S = ", no_species)
end

# randomize guild interactions

# fraction to be randomized
π = 0.15

# tally number of interactions
no_interactions = [0]
A_col_sum = sum(A,dims=2)
#no_interactions = sum(A_col_sum,dims=1)
for i = 1:no_guilds
    for j = 1:no_guilds
        if P[i,:G]>0 && P[j,:G]>0
            no_interactions[1] = no_interactions[1] + 1
        end
    end
end

println("No. of interactions = ",no_interactions[1])

# no. to be reassigned
reassign = convert(Int64,round(no_interactions[1]*π))
println("No. of interactions to be randomized = ",reassign)

# number of randomizations
for sims = 1:1
    
    # reassign in A matrix
    for i = 1:reassign
        
        # make vector of guilds and shuffle
        guild_nos = Int64[]
        for j = 1:no_guilds
            push!(guild_nos,j)
        end
        shuffle!(guild_nos)
        
        # select guild and interaction, randomize and swap one
        ones = Int64[]
        naughts = Int64[]
        random_guild = guild_nos[1]
        for k = 1:no_guilds
            if A[random_guild,k]==0
                push!(naughts,k)
            else
                if A[random_guild,k]==1
                    push!(ones,k)
                end
            end
        end
        # avoid empty guilds
        if size(ones,1)==0 || size(naughts,1)==0
            continue
        end
        # randomize and reassign
        shuffle!(ones); shuffle!(naughts)
        A[random_guild,ones[1]] = 0
        A[random_guild,naughts[1]] = 1
    end
    
    # simulate n SLNs
    s = species_level_network(A,P,γ)
    # write to file
    for j = 1:size(s,1)
        print(f1,sims)
        for k = 1:size(s,2)
            print(f1,",",s[j,k])
        end
        print(f1,"\n")
    end
    
end # end simulations

close(f1)



