### A Pluto.jl notebook ###
# v0.19.11

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
	Pkg.activate("/Users/pavanchaggar/ResearchDocs/Presentations/Roche-0622")
end

# ╔═╡ bc07c567-bf0e-449f-adcf-afb127700900
begin
	using PlutoUI
	using Plots
	using DifferentialEquations
	using HypertextLiteral: @htl
	using Connectomes
end

# ╔═╡ 91107cb3-a72e-47a7-8d21-42b2ea11e521
include("/Users/pavanchaggar/ResearchDocs/Presentations/Roche-0622/functions.jl")

# ╔═╡ c051a690-e854-436f-993f-0fa80b000a73
html"""<style>
main {
max-width: 900px;
}"""

# ╔═╡ caeee295-cecd-4261-8846-2511cbb91d41
html"<button onclick='present()'>present</button>"

# ╔═╡ 9b698c74-dbc7-4510-ae29-ead164bcf830
c = filter(Connectome("/Users/pavanchaggar/.julia/dev/Connectomes/assets/connectomes/Connectomes-hcp-scale1.xml"));

# ╔═╡ 2aae6fa9-a3f0-451a-80dd-8b119f48072d
const L = laplacian_matrix(c);

# ╔═╡ 2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
md"
# Mathematical Modelling and Inference Methods for Alzheimer's Disease

**Pavanjit Chaggar, June 2022** \
pavanjit.chaggar@maths.ox.ac.uk \
@ChaggarPavan on Twitter

DPhil student at the Mathematical Institute, University of Oxford.
Supervised by Alain Goriely, Saad Jbabdi, Stefano Magon and Gregory Klein.
"

# ╔═╡ 0f3da277-c6ca-484f-9f83-b5899a3b2d5f
md"
# Aim
#### To build up models of AD that effectively describe our data.
"

# ╔═╡ a828a333-df39-4a4a-8744-0c235fb4342e
md"
# Overview and Introduction

- Alzheimer's disease (AD)
- Mathematical models of AD
- Patient Inference Case Study
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
We want to build up models like lego. 
"

# ╔═╡ b0618ecd-e43e-4378-b90b-5f480a601749
md" 
## Structural Connectomes and Transport
"

# ╔═╡ 5c30120e-7923-4891-8f7f-b086bbf7f3e6
md"
The first important part of the modelling of $\tau$P in AD is describing **transport through the brain**. In this work, we model transport as diffusion across the structural network. We obtain the structural connectome using tractography data from HCP, processed using ProbTrackX in FSL.
"

# ╔═╡ 83539771-b2bd-4ab0-b1e5-2444323c21e9
pic("https://github.com/PavanChaggar/Presentations/blob/master/Roche-1221/assets/images/connectomes/connectome-length-free.png"; h =300, w=900)

# ╔═╡ 2b2e5e0b-7ac6-40c6-84ed-d5e12fd64e95
begin
        function NetworkDiffusion(du, u, p, t)
        du .= -p * L * u
        end
        function simulate(prob, p)
                solve(remake(prob, p=p), Tsit5())
        end;

        u1 = zeros(83)
        u1[[27, 68]] .= 0.5

        p1 = 1.0
        t_span1 = (0.0,20.0)

        prob_diffusion = ODEProblem(NetworkDiffusion, u1, t_span1, p1)

        sol_diffusion = solve(prob_diffusion, Tsit5())
end;

# ╔═╡ bf14982a-7d07-4a24-8353-bab87a06f56f
md"## Diffusion Model
"

# ╔═╡ 6df47d17-e4a5-4f8b-82d6-f8c12b8f19e9
two_cols(md"",md"
ρ = $(@bind ρ1 PlutoUI.Slider(0:0.1:5, show_value=true, default=0))")

# ╔═╡ 8309cc04-26a7-4896-bf83-76977c6dd28f
TwoColumn(
md"
Using the graph Laplacian, we can define the network heat equation, which describes diffusion across a graph. 

To the left, I have shown a simulation with an initial seeding concentration placed in the entorhinal cortex. By changing the diffusion coefficient, $\rho$, we can see how the dynamics are affected.
\
\
$\frac{dp_i}{dt} = \underbrace{\sum_{j} -\rho L_{ij} p_j}_{transport}$
"
,     
Plots.plot(simulate(prob_diffusion, ρ1), size=(450,300), labels=false, ylims=(0.0,0.5), xlims=(0.0,20.0), ylabel="Concentration")
)

# ╔═╡ 93e36f97-bd13-43c3-8d02-d756c223383c
md" ## Diffusion Model"

# ╔═╡ 3af2f496-23b4-41fc-8072-69cf83b1d2fa
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/diffusion.mp4") 

# ╔═╡ dd63aa8e-2ef9-4d18-8f2e-cda1a825efaa
begin
		prob_fkpp = ODEProblem(NetworkFKPP, u1, t_span1, [0.1,1.0])
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
The next piece we want to add to the model is **autocatalytic protein growth**, to more accurately describe toxic protein dynamics. The most simple such term we could include is a quadratic term that is bounded between $0$ and $1$. 
\
\
The effect of this quadratic term is exponential growth given a positive concentration of toxic protein that saturates as the concentration grows.
\
\
$$\frac{d p_i}{dt} = \underbrace{\sum_j -\rho L_{ij}p_j}_{transport} + \underbrace{\alpha p_i\left(1-p_i\right)}_{growth}$$
",     
Plots.plot(simulate(prob_fkpp, [ρ, α]), size=(450,300), labels=false, ylims=(0.0,1.0), xlims=(0.0,20.0), ylabel="concentration"))


# ╔═╡ 84f50b04-25d1-412e-bacf-5c0e9299eb63
md"
## FKPP Model"

# ╔═╡ 607d0291-89f3-4d4e-bb53-cc4de43de049
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/fkpp.mp4")

# ╔═╡ 696cf4fb-4687-4306-b74b-b375215d1a1f
md" 
##  Generalising the FKPP model
"

# ╔═╡ a2c1eb7c-a782-4ff3-9171-7acf1bb3d338
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/models/generalised-fkpp.png"; h = 350, w=800)

# ╔═╡ 57f7b7e2-ded0-4eac-87a4-2077b3522535
md"## Generalising the FKPP model: Local Parameters"

# ╔═╡ 15fbec7e-ae2c-4ffe-86c4-b6b1beacdfb3
two_cols(md"
Now we have transport and growth, the next piece we want to add is **regional specificity**, such as baseline SURV values and carrying capacities.

We *estimate* these using Gaussian mixture modelling of population SUVR data per region.

GMMs are fit to data from BioFinder, which has much better coverage of late stage AD subjects than ADNI, resulting in less sampling bias.",
	
pic("https://github.com/PavanChaggar/Presentations/blob/master/Roche-0622/assets/images/gmm-lEC.png"; h = 275, w=450))

# ╔═╡ 89dcf294-7b17-4aa4-8ad3-77dfd3a2d808
md" 
## Generalising the FKPP Model: Carrying Capacities
"

# ╔═╡ 724ebc25-d902-47e5-8ff4-916c77424768
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/models/carrying-capacities.png"; h = 350, w=900)

# ╔═╡ ece49802-e660-48fb-8592-f9a4098f10e8
md"
## Generalising the FKPP Model: Dynamics"

# ╔═╡ ef098338-1b67-4682-bd05-e4154e5a420f
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/exfkpp-all.mp4")

# ╔═╡ dc8da42d-afdb-423b-812e-01160ccf637a
md" 
# Summary

* Models are like lego!
* Each piece we add make the dynamics more expressive and better able to describe AD. 
"

# ╔═╡ d3a9829f-7ac4-4465-acb5-277d09cacce4
md" 
# Fitting to Subject Data
'Ok, Pavan, you've shown off your models, now tell us how they're useful?'
"

# ╔═╡ cf1e590b-e44f-4b33-bbda-cfe4ad579cb6
md" 
# Identifying Seeding Locations
* This is a **VERY** hard problem. 
* Typically not possible using the global FKPP model, since it is a model of **concentration**. It will either (1) provide trivial solutions where the region with the highest concentration will be identified as the seed. (2) Become non-identifiable due to saturation.
* The problem is more tractable using the local FKPP model, since it models SUVR with regional baseline values and carrying capacities. However, in general, seeing sites will still be non-identifiable after some nodes are saturated. 
* To do this, we perform inference using HMC and use a *horseshoe* prior to enforce sparsity. 
"

# ╔═╡ 00d6a9ac-1173-4a0d-9f3e-58e8ab6a6959
md"
## Identifying Seeding Locations
"

# ╔═╡ 3997d798-2b48-4b90-a15a-6fef3c5ecbb6
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/seeding-locations-no-divs.png"; h = 450, w=900)

# ╔═╡ be5bf2e6-96fa-4134-a6c8-3fa6a97c123c
md"
# Do Different Seeding Locations Explain Subtypes? 
"

# ╔═╡ 95d1b675-d47d-4451-a924-136157d76358
md" 
### Entorhinal Seeding
"

# ╔═╡ 77d563e6-9846-4aee-bdf5-9e01dcb6d2c6
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/localfkpp-cortical-entorhinal.mp4")

# ╔═╡ 9c8cdfbd-1848-4240-9b7d-59cc7cb9da65
md"
### Inferior Temporal Seeding
"

# ╔═╡ 7434dc23-3a6d-450f-86f1-7b5974f8801f
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/localfkpp-cortical-inferiortemporal.mp4")

# ╔═╡ 739aa309-fdec-478c-80c1-8f1efa0509bc
md" 
# Hierarchical Modelling

Using hierarchical models we can: 
* Compare parameter values across different populations, e.g. Aβ+ and Aβ- patients.
* Limit overfitting to single subjects.
* Aid identifiability, grouping information across subjects.

" 

# ╔═╡ 55d79b31-c2fe-4d5a-8776-28e7cb815666
md"
## Population Level Distributions
" 

# ╔═╡ fbfa927d-5043-4386-9589-9589f85bec1d
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/hier-dsts-wlabels.png"; h = 450, w=900)

# ╔═╡ 91f4e92c-6bba-4380-b309-f78ffd07329b
md" 
## Subject Level Distributions
" 

# ╔═╡ 6158d5cf-fe01-4c04-ad83-d921546e6f41
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/sub-dsts.png"; h = 450, w=900)

# ╔═╡ 255536ce-88eb-474b-b465-84a75edbd767
md" 
# Limitations + Conclusions 
* We can use models to describe various processes in AD. 
* ...However, the models are macro-scale descriptions of AD that don't provide a lot of mechanistic **explanation**.

* Using Bayesian inference, we can calibrate models and estimate patient disease trajectories. 
* ...However, we have very limited and noisy longitudinal patient data that makes it hard to perform inference and model comparison. 
"

# ╔═╡ 82411fe9-4773-4aea-8710-f2ae15692585
md"# Questions?"

# ╔═╡ Cell order:
# ╠═1f540848-eb08-11ec-32c6-d78736f8362e
# ╠═c051a690-e854-436f-993f-0fa80b000a73
# ╠═caeee295-cecd-4261-8846-2511cbb91d41
# ╠═bc07c567-bf0e-449f-adcf-afb127700900
# ╠═9b698c74-dbc7-4510-ae29-ead164bcf830
# ╠═2aae6fa9-a3f0-451a-80dd-8b119f48072d
# ╠═91107cb3-a72e-47a7-8d21-42b2ea11e521
# ╟─2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
# ╟─0f3da277-c6ca-484f-9f83-b5899a3b2d5f
# ╟─a828a333-df39-4a4a-8744-0c235fb4342e
# ╟─5ff7a99d-0ea0-4919-8ffc-a41ab94984fe
# ╟─c7bd3abf-e615-4acd-a0b5-24e80ecfee68
# ╟─654bdbd1-3190-45dc-9d71-a6eb9ade28c5
# ╟─1212f837-541e-48cc-9caf-3115aee37987
# ╟─b0618ecd-e43e-4378-b90b-5f480a601749
# ╟─5c30120e-7923-4891-8f7f-b086bbf7f3e6
# ╟─83539771-b2bd-4ab0-b1e5-2444323c21e9
# ╟─2b2e5e0b-7ac6-40c6-84ed-d5e12fd64e95
# ╟─bf14982a-7d07-4a24-8353-bab87a06f56f
# ╟─8309cc04-26a7-4896-bf83-76977c6dd28f
# ╟─6df47d17-e4a5-4f8b-82d6-f8c12b8f19e9
# ╟─93e36f97-bd13-43c3-8d02-d756c223383c
# ╟─3af2f496-23b4-41fc-8072-69cf83b1d2fa
# ╟─dd63aa8e-2ef9-4d18-8f2e-cda1a825efaa
# ╟─1e157396-03ae-43ff-a3a1-dec8776507e6
# ╟─ed5dfdbf-db67-47cd-8a06-dbf7c80dc336
# ╟─a11bfbd6-703b-427b-ae20-931dc40e7973
# ╟─84f50b04-25d1-412e-bacf-5c0e9299eb63
# ╟─607d0291-89f3-4d4e-bb53-cc4de43de049
# ╟─696cf4fb-4687-4306-b74b-b375215d1a1f
# ╟─a2c1eb7c-a782-4ff3-9171-7acf1bb3d338
# ╟─57f7b7e2-ded0-4eac-87a4-2077b3522535
# ╟─15fbec7e-ae2c-4ffe-86c4-b6b1beacdfb3
# ╟─89dcf294-7b17-4aa4-8ad3-77dfd3a2d808
# ╟─724ebc25-d902-47e5-8ff4-916c77424768
# ╟─ece49802-e660-48fb-8592-f9a4098f10e8
# ╟─ef098338-1b67-4682-bd05-e4154e5a420f
# ╟─dc8da42d-afdb-423b-812e-01160ccf637a
# ╟─d3a9829f-7ac4-4465-acb5-277d09cacce4
# ╟─cf1e590b-e44f-4b33-bbda-cfe4ad579cb6
# ╟─00d6a9ac-1173-4a0d-9f3e-58e8ab6a6959
# ╟─3997d798-2b48-4b90-a15a-6fef3c5ecbb6
# ╟─be5bf2e6-96fa-4134-a6c8-3fa6a97c123c
# ╟─95d1b675-d47d-4451-a924-136157d76358
# ╟─77d563e6-9846-4aee-bdf5-9e01dcb6d2c6
# ╟─9c8cdfbd-1848-4240-9b7d-59cc7cb9da65
# ╟─7434dc23-3a6d-450f-86f1-7b5974f8801f
# ╟─739aa309-fdec-478c-80c1-8f1efa0509bc
# ╟─55d79b31-c2fe-4d5a-8776-28e7cb815666
# ╟─fbfa927d-5043-4386-9589-9589f85bec1d
# ╟─91f4e92c-6bba-4380-b309-f78ffd07329b
# ╟─6158d5cf-fe01-4c04-ad83-d921546e6f41
# ╟─255536ce-88eb-474b-b465-84a75edbd767
# ╟─82411fe9-4773-4aea-8710-f2ae15692585
