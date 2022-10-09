### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

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

# ╔═╡ 2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
md"
# Mathematical Modelling for Alzheimer's Disease

**Pavanjit Chaggar, Aug 2022** \
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

# ╔═╡ b0618ecd-e43e-4378-b90b-5f480a601749
md" 
## Structural Connectomes and Transport
"

# ╔═╡ 5c30120e-7923-4891-8f7f-b086bbf7f3e6
md"
The first important part of the modelling of $\tau$P in AD is describing **transport through the brain**. In this work, we model transport as diffusion across the structural network. We obtain the structural connectome using diffusion data from HCP, processed using ProbTrackX in FSL.
"

# ╔═╡ 83539771-b2bd-4ab0-b1e5-2444323c21e9
pic("https://github.com/PavanChaggar/Presentations/blob/master/Roche-1221/assets/images/connectomes/connectome-length-free.png"; h =300, w=900)

# ╔═╡ bf14982a-7d07-4a24-8353-bab87a06f56f
md"## Diffusion Model
"

# ╔═╡ 84f50b04-25d1-412e-bacf-5c0e9299eb63
md"
## FKPP Model"

# ╔═╡ 607d0291-89f3-4d4e-bb53-cc4de43de049
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/global-fkpp.mp4")

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
* Typically not possible using the global FKPP model, since it is a model of **concentration**. It will either (1) provide trivial solutions where the region with the highest concentration will be identified as the seed or (2) Become non-identifiable due to saturation.
* The problem is more tractable using the local FKPP model, since it models SUVR with regional baseline values and carrying capacities. However, in general, seeing sites will still be non-identifiable after some nodes are saturated. 
* We perform inference using HMC and use a *horseshoe* prior to enforce sparsity. 
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
## Entorhinal Seeding
"

# ╔═╡ 77d563e6-9846-4aee-bdf5-9e01dcb6d2c6
LocalResource("/Users/pavanchaggar/Projects/model-selection/adni/visualisation/videos/localfkpp-cortical-entorhinal.mp4")

# ╔═╡ 9c8cdfbd-1848-4240-9b7d-59cc7cb9da65
md"
## Inferior Temporal Seeding
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

# ╔═╡ 54f54339-1f56-4816-9c3d-c3667aceb8d4
pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/sub-dsts.png"; h = 450, w=900)

# ╔═╡ c5825f18-6ec4-4dcb-9e9b-2ab6d9f8bd7f
md"
## Predictions: EC
"

# ╔═╡ 081aae2f-683d-40b4-beda-7079eec5cee5
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-mtlpos-ec.png"; h = 400, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-tauneg-ec.png"; h = 400, w=900)
)

# ╔═╡ ec9fa380-52c0-4c6e-82d1-88e396c4876e
md" 
## Predictions: Hippocampus
"

# ╔═╡ 0c915672-2475-4543-92a1-1220822500bf
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-mtlpos-hc.png"; h = 400, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-tauneg-hc.png"; h = 400, w=900)
)

# ╔═╡ a0fbefc4-b34f-42fc-b3aa-289888700687
md"
## Predictions: Global tau"

# ╔═╡ dd052a6e-e81a-4d1b-b7ec-0ec0c49a3168
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-globaltaupos.png"; h = 400, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pstpred-globaltauneg.png"; h = 400, w=900)
)

# ╔═╡ 781dcc71-d62d-4f52-8ab8-0926f55c9c90
md"
## Predictions: Regional change over time
"

# ╔═╡ 9597d14c-2d7d-4366-b421-ace3c1b2d082
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pred-delta-taupos.png"; h = 350, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pred-delta-tauneg.png"; h = 350, w=900)
)

# ╔═╡ 7aefd2fa-a0f3-4460-82b1-1bf5d6e3b707
md"
## Predictions: Global change over time
"

# ╔═╡ 73587dd6-0b10-44ff-8949-60eb8eea7e59
two_cols(
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pred-delta-totaltaupos.png"; h = 350, w=900),
	pic("https://github.com/PavanChaggar/pluto-presentations/blob/main/assets/images/hier-inf/pred-delta-totaltauneg.png"; h = 350, w=900)
)

# ╔═╡ 82411fe9-4773-4aea-8710-f2ae15692585
md"# Next Steps...
* Model validation on external dataset, i.e. BioFINDER? 
* Adding amyloid interaction into the model to test whether this explains increased production rates and local variations. 
* Using gene maps to add more regional information
* Coupling SUVR with atropy to model decreases in SUVR during late stage disease. 
"

# ╔═╡ Cell order:
# ╠═1f540848-eb08-11ec-32c6-d78736f8362e
# ╠═c051a690-e854-436f-993f-0fa80b000a73
# ╠═caeee295-cecd-4261-8846-2511cbb91d41
# ╠═bc07c567-bf0e-449f-adcf-afb127700900
# ╠═91107cb3-a72e-47a7-8d21-42b2ea11e521
# ╟─2e3cfdf8-2c92-422f-917e-fc5c8b2a3451
# ╟─0f3da277-c6ca-484f-9f83-b5899a3b2d5f
# ╟─a828a333-df39-4a4a-8744-0c235fb4342e
# ╟─5ff7a99d-0ea0-4919-8ffc-a41ab94984fe
# ╟─c7bd3abf-e615-4acd-a0b5-24e80ecfee68
# ╠═654bdbd1-3190-45dc-9d71-a6eb9ade28c5
# ╟─b0618ecd-e43e-4378-b90b-5f480a601749
# ╟─5c30120e-7923-4891-8f7f-b086bbf7f3e6
# ╟─83539771-b2bd-4ab0-b1e5-2444323c21e9
# ╟─bf14982a-7d07-4a24-8353-bab87a06f56f
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
# ╟─54f54339-1f56-4816-9c3d-c3667aceb8d4
# ╟─c5825f18-6ec4-4dcb-9e9b-2ab6d9f8bd7f
# ╟─081aae2f-683d-40b4-beda-7079eec5cee5
# ╟─ec9fa380-52c0-4c6e-82d1-88e396c4876e
# ╟─0c915672-2475-4543-92a1-1220822500bf
# ╟─a0fbefc4-b34f-42fc-b3aa-289888700687
# ╟─dd052a6e-e81a-4d1b-b7ec-0ec0c49a3168
# ╟─781dcc71-d62d-4f52-8ab8-0926f55c9c90
# ╟─9597d14c-2d7d-4366-b421-ace3c1b2d082
# ╟─7aefd2fa-a0f3-4460-82b1-1bf5d6e3b707
# ╟─73587dd6-0b10-44ff-8949-60eb8eea7e59
# ╟─82411fe9-4773-4aea-8710-f2ae15692585
