### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1f540848-eb08-11ec-32c6-d78736f8362e
begin
	using Pkg
	Pkg.activate("/Users/pavanchaggar/ResearchDocs/pluto-presentations")
end

# ╔═╡ bc07c567-bf0e-449f-adcf-afb127700900
begin
	using PlutoUI
	using Plots
	using DifferentialEquations
	using HypertextLiteral: @htl
	using Connectomes
	using Images
end

# ╔═╡ 91107cb3-a72e-47a7-8d21-42b2ea11e521
include("/Users/pavanchaggar/ResearchDocs/pluto-presentations/functions.jl")

# ╔═╡ c051a690-e854-436f-993f-0fa80b000a73
html"""<style>
main {
max-width: 900px;
}"""

# ╔═╡ caeee295-cecd-4261-8846-2511cbb91d41
html"<button onclick='present()'>present</button>"

# ╔═╡ 9b698c74-dbc7-4510-ae29-ead164bcf830
c = filter(Connectome(Connectomes.connectome_path()));

# ╔═╡ 2aae6fa9-a3f0-451a-80dd-8b119f48072d
const L = laplacian_matrix(c);

# ╔═╡ 2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
md"
# Mathematical Modelling and Inference Methods for Alzheimer's Disease

**Pavanjit Chaggar, Oct 2022** \
pavanjit.chaggar@maths.ox.ac.uk \
@ChaggarPavan on Twitter

DPhil student at the Mathematical Institute, University of Oxford.
Supervised by Alain Goriely, Saad Jbabdi, Stefano Magon (Roche) and Gregory Klein (Roche).
"

# ╔═╡ 2f2105f5-aa6f-42d6-be5d-d9de2efedf5d
md" 
# Aim
Can we use mathematical models of AD to describe our in-vivo human observations? 
"

# ╔═╡ 5ff7a99d-0ea0-4919-8ffc-a41ab94984fe
md"
# Alzheimer's Disease -- A Brief Summary
Alzheimer's is characterised by gradual neurodegeneration associated with the pervasive spreading of toxic proteins.

In particular, two proteins, Amyloid beta (Aβ) and tau-protein (τP), are believed to underlie and drive the development of pathology.

Here, I will focus on the role of τP, in part because it spreads very predictably and is more tightly coupled with atrophy and symptom onset.
"


# ╔═╡ c7bd3abf-e615-4acd-a0b5-24e80ecfee68
md"## Braak Stages of Tau protein
In most AD cases, τP follows a predictable pattern of spreading, starting in the entorhinal cortex before spreading through the hippocampal regions, lateral cortex and finally into the neocortex. Atrophy tends to follow the spreading pattern of $\tau$P, more so than that of Aβ."

# ╔═╡ 654bdbd1-3190-45dc-9d71-a6eb9ade28c5
pic("https://github.com/PavanChaggar/Presentations/blob/master/Roche-1221/assets/images/braak-stages.png"; h=300, w=900)

# ╔═╡ 1212f837-541e-48cc-9caf-3115aee37987
md" 
# Modelling on Brain Networks! 
"

# ╔═╡ b0618ecd-e43e-4378-b90b-5f480a601749
md" 
## Structural Connectomes and Transport
"

# ╔═╡ 5c30120e-7923-4891-8f7f-b086bbf7f3e6
md"
The first important part of the modelling of $\tau$P in AD is describing **transport through the brain**. In this work, we model transport as diffusion across the structural network. We obtain the structural connectome using diffusion data from HCP, processed using ProbTrackX in FSL.
"

# ╔═╡ 83539771-b2bd-4ab0-b1e5-2444323c21e9
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/connectomes/weightednetwork.png"; h =300, w=900)

# ╔═╡ 6bd79980-f152-4dcf-b3c9-4a53f4578674
md" 
## Diffusion Model 
" 

# ╔═╡ c85decbf-d539-4bf9-aa6e-c61bd8bbeceb
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/diffusion-cortical.mp4")

# ╔═╡ 57a48ec5-f690-467b-86f1-e22b0fb4e73c
md"
## Protein Growth
"

# ╔═╡ 1525f746-e2e3-4776-b458-22462270fbda
two_cols(pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/models/heterodimerkinetics.png"; h=200, w=400),
md"
```math
\begin{align}
    \dfrac{\mathrm{d}p_i}{\mathrm{d}t} &= -\rho L_{ij} p_j +  k_0 - k_1 p_i - k_{12}p_i\tilde{p}_i \\
    \dfrac{\mathrm{d}\tilde{p}_{i}}{\mathrm{d}t} &= -\rho L_{ij} p_j \phantom{+ k_0,} - \tilde{k}_1 \tilde{p}_i + k_{12}p_i\tilde{p}_i
\end{align}
```
"
)

# ╔═╡ dd63aa8e-2ef9-4d18-8f2e-cda1a825efaa
begin
		function fkpp(du, u, p, t; L = L)
        	du .= -p[1] * L * u .+ p[2] .* u .* (1 .- u)
		end
	  	function simulate(prob, p)
                solve(remake(prob, p=p), Tsit5())
        end;

 		u1 = zeros(83)
        u1[[27, 68]] .= 0.5

        p1 = 1.0
        t_span1 = (0.0,20.0)
		prob_fkpp = ODEProblem(fkpp, u1, t_span1, [0.1,1.0])
        sol_fkpp = solve(prob_fkpp, Tsit5())
end;

# ╔═╡ 1e157396-03ae-43ff-a3a1-dec8776507e6
md" 
## FKPP Model
"

# ╔═╡ a11bfbd6-703b-427b-ae20-931dc40e7973
two_cols(md"",
md"
ρ = $(@bind ρ Slider(0:0.1:5, show_value=true, default=0)) \
α = $(@bind α Slider(-3:0.1:3, show_value=true, default=0))
")

# ╔═╡ ed5dfdbf-db67-47cd-8a06-dbf7c80dc336
TwoColumn(
md"
The FKPP model is a reaction-diffusion model, derived from a kinetic model of the toxic prion process.

The model describes protein **transport** across a network and local concentration **growth**. 
\
\
The second term describes exponential growth given a positive concentration of toxic protein that saturates toward $1.0$ as the concentration grows.
\
\
$$\frac{d p_i}{dt} = \underbrace{\sum_j -\rho L_{ij}p_j}_{transport} + \underbrace{\alpha p_i\left(1-p_i\right)}_{growth}$$
",     
Plots.plot(simulate(prob_fkpp, [ρ, α]), size=(425,250), labels=false, ylims=(0.0,1.0), xlims=(0.0,20.0), ylabel="concentration", linewidth=2))


# ╔═╡ 84f50b04-25d1-412e-bacf-5c0e9299eb63
md"
## FKPP Model"

# ╔═╡ 607d0291-89f3-4d4e-bb53-cc4de43de049
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/global-fkpp.mp4")

# ╔═╡ 679afb1d-8d72-41b1-9648-2ed1bc2c5027
md" 
# Incorporating Regional Information 

Can we incorporate regional information about $\tau P$ SUVR dynamics?
" 

# ╔═╡ e69b89f8-c96c-48af-b1b3-81f2e9ca072c
md" 
## Describing Heterogeneity in $\tau P$ SUVR
" 

# ╔═╡ f18e7628-f03a-4e46-a086-b68b2422d400
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/braak/scholl-pet.png"; h=350, w=750)

# ╔═╡ 4062217d-b34e-4e3b-ab5c-ef2fc9de34e2
cite("Schöll, M., Lockhart, S.N., Schonhaut, D.R., O’Neil, J.P., Janabi, M., Ossenkoppele, R., Baker, S.L., Vogel, J.W., Faria, J., Schwimmer, H.D. and Rabinovici, G.D., 2016. PET imaging of tau deposition in the aging human brain. Neuron, 89(5), pp.971-982")

# ╔═╡ 696cf4fb-4687-4306-b74b-b375215d1a1f
md" 
##  Generalising the FKPP model
"

# ╔═╡ c2c35600-219c-4593-a8bf-6c296ac1bda4
plot(load("/Users/pavanchaggar/ResearchDocs/pluto-presentations/assets/images/models/models.png"), showaxis=false, ticks=false, size=(750, 350))

# ╔═╡ 57f7b7e2-ded0-4eac-87a4-2077b3522535
md"## Generalising the FKPP model: Local Parameters"

# ╔═╡ 15fbec7e-ae2c-4ffe-86c4-b6b1beacdfb3
two_cols(md"
Now we have transport and growth, the next piece we want to add is **regional specificity**, such as baseline SURV values and carrying capacities.

We *estimate* these using Gaussian mixture modelling of population SUVR data per region.

GMMs are fit to data from BioFinder, which has much better coverage of late stage AD subjects than ADNI, resulting in less sampling bias.",
	
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/gmm-lIT.png"; h = 250, w=450))

# ╔═╡ ece49802-e660-48fb-8592-f9a4098f10e8
md"
## Generalising the FKPP Model: Dynamics"

# ╔═╡ ef098338-1b67-4682-bd05-e4154e5a420f
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/localfkpp-cortical-entorhinal.mp4")

# ╔═╡ 739aa309-fdec-478c-80c1-8f1efa0509bc
md" 
# Hierarchical Modelling

Using hierarchical models we can: 
* Compare parameter values across different populations. Here we compare Aβ-, Aβ+τP- and Aβ+τP+ patients.
* Limit overfitting to single subjects.
* Aid identifiability, grouping information across subjects.

" 

# ╔═╡ 55d79b31-c2fe-4d5a-8776-28e7cb815666
md"
## Population Level Distributions
" 

# ╔═╡ fbfa927d-5043-4386-9589-9589f85bec1d
plot(load("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/hier-inf/hier-dsts.pdf"), showaxis=false, ticks=false, size=(1000, 500))

# ╔═╡ 91f4e92c-6bba-4380-b309-f78ffd07329b
md" 
## Subject Level Distributions
" 

# ╔═╡ 54f54339-1f56-4816-9c3d-c3667aceb8d4
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/sub-dsts.png"; h =350, w=750)

# ╔═╡ 22931974-4961-4558-89fb-316f2c207484
md"
## Predictions
"

# ╔═╡ d841bad5-3244-4634-910c-2a4bf8c83dec
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pred-taupos.png"; h=350, w=700, hspace=50)

# ╔═╡ ec9fa380-52c0-4c6e-82d1-88e396c4876e
md" 
## Predictions: Hippocampus
"

# ╔═╡ 0c915672-2475-4543-92a1-1220822500bf
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-taupos-Right-Hippocampus.png"; h =400, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-tauneg-Right-Hippocampus.png"; h =330, w=900)
)

# ╔═╡ 82411fe9-4773-4aea-8710-f2ae15692585
md"# Next Steps...
* Model validation on external dataset, i.e. BioFINDER? 
* Adding amyloid interaction into the model to test whether this explains increased production rates and local variations. 
* Using gene maps to add more regional information.
* Coupling SUVR with atropy to model decreases in SUVR during late stage disease. 
"

# ╔═╡ Cell order:
# ╠═1f540848-eb08-11ec-32c6-d78736f8362e
# ╠═c051a690-e854-436f-993f-0fa80b000a73
# ╠═caeee295-cecd-4261-8846-2511cbb91d41
# ╠═bc07c567-bf0e-449f-adcf-afb127700900
# ╠═91107cb3-a72e-47a7-8d21-42b2ea11e521
# ╠═9b698c74-dbc7-4510-ae29-ead164bcf830
# ╠═2aae6fa9-a3f0-451a-80dd-8b119f48072d
# ╟─2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
# ╟─2f2105f5-aa6f-42d6-be5d-d9de2efedf5d
# ╟─5ff7a99d-0ea0-4919-8ffc-a41ab94984fe
# ╟─c7bd3abf-e615-4acd-a0b5-24e80ecfee68
# ╟─654bdbd1-3190-45dc-9d71-a6eb9ade28c5
# ╟─1212f837-541e-48cc-9caf-3115aee37987
# ╟─b0618ecd-e43e-4378-b90b-5f480a601749
# ╟─5c30120e-7923-4891-8f7f-b086bbf7f3e6
# ╟─83539771-b2bd-4ab0-b1e5-2444323c21e9
# ╟─6bd79980-f152-4dcf-b3c9-4a53f4578674
# ╟─c85decbf-d539-4bf9-aa6e-c61bd8bbeceb
# ╟─57a48ec5-f690-467b-86f1-e22b0fb4e73c
# ╟─1525f746-e2e3-4776-b458-22462270fbda
# ╟─dd63aa8e-2ef9-4d18-8f2e-cda1a825efaa
# ╟─1e157396-03ae-43ff-a3a1-dec8776507e6
# ╟─ed5dfdbf-db67-47cd-8a06-dbf7c80dc336
# ╟─a11bfbd6-703b-427b-ae20-931dc40e7973
# ╟─84f50b04-25d1-412e-bacf-5c0e9299eb63
# ╟─607d0291-89f3-4d4e-bb53-cc4de43de049
# ╟─679afb1d-8d72-41b1-9648-2ed1bc2c5027
# ╟─e69b89f8-c96c-48af-b1b3-81f2e9ca072c
# ╟─f18e7628-f03a-4e46-a086-b68b2422d400
# ╟─4062217d-b34e-4e3b-ab5c-ef2fc9de34e2
# ╟─696cf4fb-4687-4306-b74b-b375215d1a1f
# ╟─c2c35600-219c-4593-a8bf-6c296ac1bda4
# ╟─57f7b7e2-ded0-4eac-87a4-2077b3522535
# ╟─15fbec7e-ae2c-4ffe-86c4-b6b1beacdfb3
# ╟─ece49802-e660-48fb-8592-f9a4098f10e8
# ╟─ef098338-1b67-4682-bd05-e4154e5a420f
# ╟─739aa309-fdec-478c-80c1-8f1efa0509bc
# ╟─55d79b31-c2fe-4d5a-8776-28e7cb815666
# ╟─fbfa927d-5043-4386-9589-9589f85bec1d
# ╟─91f4e92c-6bba-4380-b309-f78ffd07329b
# ╟─54f54339-1f56-4816-9c3d-c3667aceb8d4
# ╟─22931974-4961-4558-89fb-316f2c207484
# ╟─d841bad5-3244-4634-910c-2a4bf8c83dec
# ╟─ec9fa380-52c0-4c6e-82d1-88e396c4876e
# ╟─0c915672-2475-4543-92a1-1220822500bf
# ╟─82411fe9-4773-4aea-8710-f2ae15692585
