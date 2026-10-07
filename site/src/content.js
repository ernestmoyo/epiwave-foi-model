// Reading, derivations, worked examples and questions for each chapter.
// Numbers quoted here come from R/epiwave-foi-model.R and the 50-replicate study.
window.EPIWAVE_CONTENT = {
  chapters: [
    {
      id: "design",
      short: "The design",
      title: "Mechanism first, statistics second",
      lede: "Why the model has two stages, and what the offset buys you.",
      minutes: 8,
      reading: [
        "A standard geostatistical map of malaria asks one question of the data: where is incidence high? It answers with an intercept and a spatial random field, and it has no idea why. Mosquitoes, bednets and temperature never enter the model; any pattern they create has to be rediscovered from the case counts.",
        "EpiWave-FOI starts from the opposite end. A Ross-Macdonald model, driven by fixed entomological inputs, predicts the infection incidence each site should have. That prediction, I*, goes into the statistical model as a fixed offset on the log scale. The Gaussian process then only has to explain where the mechanism is wrong.",
        "The split into two stages is what keeps this tractable. The differential equations are solved once per site, before any sampling, with parameters held fixed. Inference never touches the ODE, so fitting costs the same as an ordinary geostatistical model.",
        "Because the prior on the residual field is centred on zero, the model defaults to the mechanistic prediction: with no data, I = e<sup>α</sup> I*. Data pull the map away from I* only where they disagree with it. Setting I* = 0 (dropping the offset) recovers the standard geostatistical model, which is the comparison the simulation study makes."
      ],
      derivations: [
        {
          title: "Why the model defaults to the mechanism",
          steps: [
            "Write the latent incidence as log I<sub>s,t</sub> = α + log I*<sub>s,t</sub> + ε<sub>s,t</sub>.",
            "Exponentiate: I<sub>s,t</sub> = e<sup>α</sup> · I*<sub>s,t</sub> · e<sup>ε<sub>s,t</sub></sup>.",
            "The GP prior has mean zero, so before seeing data the most likely field is ε = 0 everywhere.",
            "Then I<sub>s,t</sub> = e<sup>α</sup> I*<sub>s,t</sub>: incidence is proportional to the mechanistic prediction, and α sets the constant of proportionality.",
            "Any departure from proportionality must be paid for in the GP prior, so the model only uses ε where the data demand it."
          ]
        },
        {
          title: "What I* = 0 means",
          steps: [
            "\"Setting I* = 0\" means removing the offset term, not putting log 0 into the model.",
            "The model becomes log I<sub>s,t</sub> = α + ε<sub>s,t</sub>: an intercept plus a space-time field.",
            "This is the standard geostatistical model. Everything else (priors, likelihoods, sampler) stays identical, so any difference between the two fits is due to the mechanistic information alone."
          ]
        }
      ],
      examples: [
        {
          title: "Reading the offset at one site",
          body: "Replicate 1 of the study (run with earlier, higher illustrative values), site 1, month 24: the ODE gives I* ≈ 0.093 infections per person per day. If α = 0.11 and the fitted residual is ε = −0.4, then I = e<sup>0.11</sup> × 0.093 × e<sup>−0.4</sup> ≈ 1.12 × 0.093 × 0.67 ≈ 0.070. The mechanism over-predicts this site-month by about a third, and the GP absorbs the difference."
        }
      ],
      selftest: [
        { q: "Why is the ODE solved only once, rather than inside the MCMC?", a: "Its parameters are fixed inputs, not unknowns. Nothing the sampler changes affects the ODE solution, so solving it again would return the same I*.", depth: "core" },
        { q: "What does the GP represent in this model, in one sentence?", a: "The space-time pattern in log incidence that the mechanistic prediction gets wrong.", depth: "core" },
        { q: "If the GP variance were forced to zero, what would the model assume?", a: "That incidence is exactly proportional to I* everywhere, i.e. that the mechanism explains all spatial and temporal variation up to a constant.", depth: "deeper" }
      ]
    },
    {
      id: "stage1",
      short: "Stage 1",
      title: "Ross-Macdonald and I*",
      lede: "From mosquito numbers to a predicted infection rate.",
      minutes: 10,
      reading: [
        "Stage 1 follows two proportions over time at each site: x, the share of people infected, and z, the share of mosquitoes infectious. Mosquitoes bite at rate a; there are m mosquitoes per person; a bite from an infectious mosquito infects a person with probability b, and a bite on an infected person infects a mosquito with probability c. People recover at rate r; mosquitoes die at rate g.",
        "Only m, a and g vary in space and time. b, c and r are held fixed everywhere because they are not expected to vary much between places. m varies through the year with two rainy seasons; a and g are constant stand-ins until mapped surfaces replace them.",
        "The quantity Stage 2 uses is the force of infection: the rate at which a susceptible person acquires infection, I* = m · a · b · z, per day. It is a rate, not a count; population enters later, in the case likelihood. Incidence among the whole population is I* (1 − x). The model treats I* as proportional to incidence, so the (1 − x) factor is absorbed by α and ε. In the simulation x is at most about 0.41, so the two differ by up to 1.7 times.",
        "The fixed values come from the literature: b = c = 0.5 (Smith & McKenzie 2004; Smith et al. 2007), infections that last about six months, r = 1/180 per day (Smith et al. 2005), and a 10-day extrinsic incubation period (Gething et al. 2011; Stopard et al. 2021). With m = 0.2, a = 0.3 and g = 0.1, R₀ is about 3.3 and prevalence before nets is 0.48, close to Mozambique's 2018 survey (39% nationally). I* is then about one infection per person per year. Vector Atlas surfaces will replace m and a."
      ],
      derivations: [
        {
          title: "The two equations",
          steps: [
            "People: infectious bites arrive at m a z per person per day; a fraction b infects. Only uninfected people (1 − x) can be infected; infected people recover at rate r.",
            "dx/dt = m a b z (1 − x) − r x",
            "Mosquitoes: each bites people at rate a; a fraction x of people are infected and a fraction c of those bites infect. Only uninfected mosquitoes (1 − z) can be infected; they die at rate g.",
            "dz/dt = a c x (1 − z) − g z"
          ]
        },
        {
          title: "R₀ for this Ross-Macdonald model",
          steps: [
            "One infected person, in a fully susceptible population, infects mosquitoes at rate m a c (m mosquitoes per person, each biting at rate a, infected with probability c), for an average of 1/r days.",
            "Mosquitoes infected per person: m a c / r.",
            "Each infected mosquito infects people at rate a b for an average of 1/g days: a b / g people.",
            "R₀ = (m a c / r) × (a b / g) = m a² b c / (g r).",
            "Without an extrinsic incubation period this is the whole story. With one of n days, only mosquitoes that survive it can transmit. The model splits it into four exposed stages, so R₀ is multiplied by (1 + g n / 4)<sup>−4</sup>, close to the exact e<sup>−g n</sup>. At g = 0.1 and n = 10 that is 0.41 (exact 0.37), and I* falls by more than R₀ does."
          ]
        },
        {
          title: "Equilibrium prevalence",
          steps: [
            "Set both derivatives to zero. From the mosquito equation: z* = a c x / (a c x + g).",
            "Put z* into the human equation and solve for x: r (a c x + g) = m a b · a c (1 − x).",
            "Collect terms: x* = (m a² b c − r g) / (a c (m a b + r)).",
            "Divide top and bottom by r g: x* = (R₀ − 1) / (R₀ + a c / g).",
            "x* > 0 only when R₀ > 1: the infection persists only if each case replaces itself."
          ]
        }
      ],
      examples: [
        {
          title: "The simulation's baseline values",
          body: "m = 0.2, a = 0.3, b = c = 0.5, g = 0.1, r = 1/180, EIP 10 days. Surviving the EIP: (1 + 0.1 × 10 / 4)<sup>−4</sup> ≈ 0.41. R₀ = 0.2 × 0.09 × 0.25 × 0.41 / (0.1 / 180) ≈ 3.3. Solving the equilibrium gives x* ≈ 0.48 and z* ≈ 0.17, so I* = 0.2 × 0.3 × 0.5 × 0.17 ≈ 0.005 per person per day, about 1.9 infections per person per year. A homogeneous Ross-Macdonald model reaches this prevalence with an EIR of only about 4 per year, below field estimates; Smith et al. (2005) attribute the gap to heterogeneous biting."
        },
        {
          title: "Every simulated site has the same I*",
          body: "In simulate_epiwave_data() each site gets the same seasonal m, the same a and g, and the same ITN scale-up. The ODE therefore returns the same I* at all ten sites: the spread across sites in any month is exactly 0. The mechanistic prediction carries seasonal information but no spatial information."
        }
      ],
      selftest: [
        { q: "Why is I* a rate and not a count?", a: "It is infectious bites that infect, per person per day (m a b z). Converting to a count needs population, which enters the Poisson case likelihood separately.", depth: "core" },
        { q: "Which parameter appears squared in R₀, and why?", a: "a, the biting rate. Transmission needs two bites: one to infect the mosquito, one for the mosquito to infect a person.", depth: "core" },
        { q: "If g doubles, what happens to R₀?", a: "It halves (R₀ ∝ 1/g). With an extrinsic incubation period it falls faster, because e<sup>−gn</sup> also shrinks.", depth: "deeper" }
      ]
    },
    {
      id: "interventions",
      short: "Interventions",
      title: "Interventions act on the entomology",
      lede: "Bednets and spraying change m, a and g before the ODE runs.",
      minutes: 6,
      reading: [
        "Interventions enter the model upstream: they change the entomological parameters, and the ODE turns those changes into a change in I*. The statistical model is not told about interventions directly.",
        "Each effect is a multiplier that equals 1 at zero coverage. An ITN kills some mosquitoes that try to feed (lower m), stops some feeding (lower a) and raises mosquito mortality (higher g). IRS raises mortality and slightly reduces feeding. Insecticide resistance lowers the effective ITN coverage by the factor 1 − 0.46 (1 − q), where q is the fraction susceptible in bioassays: the relationship used in the Vector Atlas insecticide-resistance cube. A fully resistant population keeps 54% of the net effect, because nets still block bites. IRS has its own susceptibility, since IRS insecticides are mostly not pyrethroids.",
        "The effect sizes are illustrative and follow Griffin et al. (2010) and Bhatt et al. (2015). The 0.46 was estimated on the prevalence scale; a hut-trial mapping from bioassays to killing and deterrence is the planned replacement. The current version has no human-behaviour component: how and when people are bitten, and how that changes with nets, is the next piece of Stage 1 to build."
      ],
      derivations: [
        {
          title: "The ITN multipliers",
          steps: [
            "Let n be ITN coverage and q the fraction of vectors susceptible in bioassays.",
            "Effective coverage: n′ = n · (1 − 0.46 (1 − q)).",
            "m → m · (1 − n′ k), where k is the kill rate among net encounters.",
            "a → a · (1 − n′ f), where f is the feeding inhibition.",
            "g → g · (1 + n′ μ), where μ is the mortality boost.",
            "At n = 0 every multiplier is 1. At full resistance (q = 0) nets keep 54% of their effect, because they still block bites."
          ]
        }
      ],
      examples: [
        {
          title: "70% ITN coverage, 80% susceptible",
          body: "Effective coverage n′ = 0.7 × (1 − 0.46 × 0.2) = 0.636. m is multiplied by 1 − 0.636 × 0.5 = 0.682; a by 1 − 0.636 × 0.3 = 0.809; g by 1 + 0.636 × 0.3 = 1.191. Higher mortality also means fewer mosquitoes survive the 10-day EIP (0.35 instead of 0.41). R₀ scales by 0.682 × 0.809² × 0.86 / 1.191 ≈ 0.32, so the baseline 3.3 becomes about 1.07: barely above 1, which is why prevalence keeps falling through the simulation."
        }
      ],
      selftest: [
        { q: "Why put interventions into m, a and g rather than into the statistical model?", a: "So the mechanism, not a regression coefficient, carries their effect, and the same assumptions used in transmission models apply. The GP then captures what the mechanism gets wrong, including wrong intervention effects.", depth: "core" },
        { q: "Which ITN effect reduces R₀ most, per unit change?", a: "The feeding inhibition on a, because a enters R₀ squared.", depth: "deeper" }
      ]
    },
    {
      id: "stage2",
      short: "Stage 2",
      title: "The residual field: GP in space, AR(1) in time",
      lede: "How the model describes where the mechanism is wrong.",
      minutes: 9,
      reading: [
        "The residual ε is built in two steps. Each month, a fresh spatial pattern f<sub>t</sub> is drawn from a Gaussian process with a Matérn 5/2 kernel: nearby sites get similar values, with φ setting how far the similarity reaches. Then the months are chained together: ε<sub>t</sub> = θ ε<sub>t−1</sub> + f<sub>t</sub>. With θ near 1 a site that is over-predicted this month stays over-predicted next month.",
        "This separable structure keeps the latent space small: one spatial draw of size n<sub>sites</sub> per month. A single GP over (longitude, latitude, time) would create n<sub>sites</sub> × n<sub>times</sub> latent variables and the sampler could not move through it.",
        "AR(1) is applied as a single matrix multiply rather than a loop, which TensorFlow runs quickly. That code is taken unchanged from epiwave.mapping.",
        "One consequence matters for the priors. The data inform the overall variability of ε, its stationary variance τ² = σ²/(1 − θ²), much more strongly than they inform θ. The posterior is therefore a long ridge along σ² = τ²(1 − θ²), which the sampler crosses slowly. The model now samples τ² directly and derives σ² from it."
      ],
      derivations: [
        {
          title: "Stationary variance of the AR(1) field",
          steps: [
            "ε<sub>t</sub> = θ ε<sub>t−1</sub> + f<sub>t</sub>, with Var(f<sub>t</sub>) = σ² and f<sub>t</sub> independent of the past.",
            "Var(ε<sub>t</sub>) = θ² Var(ε<sub>t−1</sub>) + σ².",
            "At stationarity Var(ε<sub>t</sub>) = Var(ε<sub>t−1</sub>) = τ².",
            "τ² = θ² τ² + σ², so τ² = σ² / (1 − θ²).",
            "Note: the model starts at ε<sub>1</sub> = f<sub>1</sub>, with variance σ², and approaches τ² over the first few months."
          ]
        },
        {
          title: "Why (σ², θ) form a ridge",
          steps: [
            "If the data inform τ² much more strongly than θ, the likelihood changes only slowly along the curve σ² = τ² (1 − θ²).",
            "Moving along that curve trades a larger innovation variance for a smaller correlation, or the reverse, without changing what the data see.",
            "A flat prior on θ adds nothing to break the tie, so chains settle at different points along it: R-hat of 2 to 4 for σ² and θ in the 50-replicate study.",
            "Sampling τ² directly puts the data's information on one parameter. A Beta(2, 2) prior on θ keeps it away from 0 and 1."
          ]
        },
        {
          title: "AR(1) as one matrix multiply",
          steps: [
            "Unroll the recursion: ε<sub>t</sub> = Σ<sub>i=0</sub><sup>t−1</sup> θ<sup>i</sup> f<sub>t−i</sub>.",
            "Build a lower-triangular matrix R with R<sub>t,s</sub> = θ<sup>t−s</sup> for s ≤ t and 0 otherwise.",
            "Then ε = R f for all months at once: one matrix product instead of a loop."
          ]
        }
      ],
      examples: [
        {
          title: "The study's true field",
          body: "σ = 0.6 and θ = 0.75 give σ² = 0.36 and a stationary variance τ² = 0.36 / (1 − 0.5625) ≈ 0.82. A site one standard deviation above zero has incidence e<sup>0.91</sup> ≈ 2.5 times its mechanistic prediction."
        },
        {
          title: "How far spatial correlation reaches",
          body: "With φ = 3 on coordinates scaled to the unit square, the Matérn 5/2 correlation between the two farthest sites of replicate 1 (distance 1.36) is 0.86. Every pair of sites has correlation at least 0.86, so the data can barely pin φ down. That is why φ is not recovered in the study."
        }
      ],
      selftest: [
        { q: "What does θ = 0 mean for the residual field?", a: "Each month's residual pattern is independent of the last: no persistence in time.", depth: "core" },
        { q: "Why not fit one GP over space and time together?", a: "It creates one latent variable per site-month, which HMC could not explore in practice. The separable GP + AR(1) keeps the latent dimension at the number of sites per month.", depth: "core" },
        { q: "If τ² is well identified but θ is not, what will the trace plots of σ² and θ look like?", a: "They will move together along a curve, slowly, and chains started in different places will disagree: the signature of a ridge.", depth: "deeper" }
      ]
    },
    {
      id: "likelihood",
      short: "Likelihoods",
      title: "Two data streams, two unknowns",
      lede: "Why prevalence must depend on I, and what that buys.",
      minutes: 9,
      reading: [
        "Case counts are modelled as Poisson(γ · I · N): γ is the reporting rate and N the population. Cases alone see only the product e<sup>α</sup> · γ. Doubling e<sup>α</sup> and halving γ predicts exactly the same counts, so the two parameters form a ridge.",
        "Prevalence surveys break the tie. The share of people testing positive depends on how many were infected recently, through I, and not on reporting. Surveys therefore inform α on its own, and cases then inform γ.",
        "This only works if survey prevalence is computed from I, the GP-adjusted incidence. An earlier version used the ODE's prevalence x, which contains neither α nor ε; the survey data then said nothing about α and the ridge stayed.",
        "Prevalence is built by summing recent infections, weighted by the chance a person infected d days ago still tests positive, q(d). epiwave.mapping uses that sum directly. This model uses p = 1 − exp(−Σ I q), the chance of at least one detectable infection. The two agree when the sum is small, and the bounded form cannot exceed 1. At the simulation's incidence they differ by about 3%. One open question: with infections that last six months, a 30-day detectability window gives survey prevalence of about 5% where the ODE says about 40%. Whether q should follow the infection's duration is for the next check-in."
      ],
      derivations: [
        {
          title: "Why cases alone cannot separate α and γ",
          steps: [
            "Expected cases: μ = γ · I · N = γ · e<sup>α</sup> · I* · e<sup>ε</sup> · N.",
            "Only the product γ e<sup>α</sup> appears. Replace (α, γ) by (α + c, γ e<sup>−c</sup>) for any c: μ is unchanged.",
            "So the case likelihood is flat along the curve γ e<sup>α</sup> = constant. The posterior is a ridge unless something else informs α."
          ]
        },
        {
          title: "From incidence to prevalence",
          steps: [
            "A person tested on day t is positive if they were infected on some earlier day t − d and are still detectable: probability q(d).",
            "Expected detectable infections per person: Λ<sub>t</sub> = Σ<sub>d</sub> I<sub>t−d</sub> q(d). On a monthly grid this becomes Λ<sub>t</sub> = Σ<sub>k</sub> I<sub>t−k</sub> w<sub>k</sub>, with w<sub>k</sub> the daily kernel integrated over each month.",
            "Linear form (epiwave.mapping): p<sub>t</sub> = Λ<sub>t</sub>.",
            "Bounded form (this model): treat detectable infections as Poisson with mean Λ<sub>t</sub>; the chance of at least one is p<sub>t</sub> = 1 − e<sup>−Λ<sub>t</sub></sup>.",
            "Expand: 1 − e<sup>−Λ</sup> = Λ − Λ²/2 + … The forms agree to first order; the bounded one stays below 1."
          ]
        },
        {
          title: "Why α in the I* = 0 model is a different quantity",
          steps: [
            "With the offset: log I = α + log I* + ε.",
            "Without it, the same incidence must be written log I = α′ + ε′.",
            "Averaging over sites and months: α′ ≈ α + mean(log I*), with the rest of log I* pushed into ε′.",
            "In replicate 1, mean(log I*) = −1.92. The I* = 0 model's α is therefore expected to sit about 1.92 below the true α. Its \"bias\" of −1.92 is this identity, not a failure.",
            "Likewise ε′ must absorb the variation in log I*, so its variance is inflated. Comparing α or σ² across the two models compares different quantities."
          ]
        }
      ],
      examples: [
        {
          title: "Detectability weights",
          body: "q(d) rises over the first few days and decays to 0 by day 30. Integrated over 30-day steps the weights are w<sub>0</sub> ≈ 10.7 and w<sub>1</sub> ≈ 6.8: a person infected this month counts for about 10.7 days of detectability, one infected last month for 6.8."
        },
        {
          title: "When the two forms part",
          body: "At the simulation's average incidence, I = 0.003 per person per day, Λ = 0.003 × (10.7 + 6.8) ≈ 0.05: linear 0.053, bounded 0.051, a 3% difference. They only part at high incidence: at I = 0.15 per day Λ ≈ 2.6, the linear form gives an impossible 2.6 and the bounded form 0.93."
        }
      ],
      selftest: [
        { q: "What would happen to the α–γ posterior if every survey were removed?", a: "It would collapse onto a ridge where γ e<sup>α</sup> is constant; only their product would be identified.", depth: "core" },
        { q: "Why did using the ODE's x in the survey likelihood fail to identify α?", a: "x does not contain α or ε, so the survey likelihood did not change when α changed. It added no information about α.", depth: "core" },
        { q: "When would the linear and bounded prevalence forms give noticeably different fits?", a: "When Λ is not small, i.e. at high incidence or with long detectability. At Λ ≈ 0.1 they differ by about 5%, and by less below that.", depth: "deeper" }
      ]
    },
    {
      id: "study",
      short: "The study",
      title: "What the simulation study shows",
      lede: "Fifty replicates, two models, one honest reading.",
      minutes: 8,
      reading: [
        "The study simulates fifty datasets from the full model (10 sites, 48 monthly steps, surveys at 30% of site-months) and fits each twice: with the I* offset and with I* = 0. Each fit used 4 chains of 2,000 samples after 2,000 warmup. That study used earlier illustrative values (7-day infections, b = c = 0.8, no EIP); the simulator now uses literature values, so the study needs rerunning.",
        "The headline comparison is the map. Median map error (RMSE of the latent incidence) is 2.4% lower with I*, and the I* model wins in 41 of the 50 replicates. That is a real but small gain.",
        "Three features of the design explain why it is small. I* is identical at every site, so it carries no spatial information. The map is scored only where data exist, and every site has cases every month, which a GP alone can fit. And surveys at 30% of site-months are far more frequent than real programmes, where a national survey comes roughly every three years.",
        "Convergence is the other caveat. α and γ mix acceptably (median R-hat 1.12 or below). φ, σ² and θ do not, in either model: median R-hat between 2 and 4. Until the sampler converges, those parameters say nothing about the model. The priors were reparameterised after this study and have not yet been tested at this scale.",
        "The redesign follows from this: different entomology at each site, the map scored at held-out sites, realistic survey frequency, and a sparse-survey scenario."
      ],
      derivations: [
        {
          title: "Why held-out sites are the right test",
          steps: [
            "At an observed site, the GP can bend the fit to the data whatever the offset, so both models fit well.",
            "At a held-out site, the I* = 0 model can only borrow from neighbours through the spatial correlation.",
            "The I* model also has the mechanistic prediction at that site. If I* carries spatial signal, that is information the I* = 0 model does not have.",
            "So the benefit of vector information should appear where data are missing, which is exactly what a map is for."
          ]
        }
      ],
      examples: [
        {
          title: "Replicate 1",
          body: "In replicate 1 the two models are almost indistinguishable: map RMSE 0.0151 with I* and 0.0150 with I* = 0. The I* = 0 model is marginally better here. Over all fifty replicates the I* model wins 41 times."
        }
      ],
      selftest: [
        { q: "Why is \"α coverage 76% vs 0%\" not evidence for the offset?", a: "Without the offset, α also absorbs mean(log I*) = −1.92, so it estimates a different quantity. The 0% coverage is guaranteed by construction.", depth: "core" },
        { q: "If I* is the same at every site, what can the offset still help with?", a: "Seasonality and the intervention trend over time: the shared temporal shape of incidence.", depth: "core" },
        { q: "What must be true before the GP hyperparameter results can be interpreted?", a: "The chains must converge (R-hat near 1, adequate ESS) for those parameters.", depth: "deeper" }
      ]
    }
  ],

  // Open questions for the next check-in. Ticked state is kept in this browser only.
  ask: [
    { id: "q-link", ch: "likelihood", q: "Prevalence link: keep the bounded 1 − exp(−Σ I·q), or use the linear Σ I·q as in epiwave.mapping?" },
    { id: "q-headline", ch: "study", q: "Is map accuracy at held-out sites the right headline metric for the with/without I* comparison?" },
    { id: "q-eip", ch: "stage1", q: "EIP and temperature-dependent c: build my own, or reuse VCOM or Stopard–Churcher code?" },
    { id: "q-sparse", ch: "study", q: "If surveys become rare, is prevalence-only fitting the next scenario to simulate?" },
    { id: "q-behaviour", ch: "interventions", q: "For biting rate × human behaviour × interventions, which data on where and when people are bitten should anchor the first version?" },
    { id: "q-paper", ch: "study", q: "For the paper: rerun the 50 replicates under the current priors, or report the earlier run with its own prior table?" },
    { id: "q-foi", ch: "stage1", q: "I* = m·a·b·z is the force of infection; incidence is I*(1 − x). Is treating I* as proportional to incidence acceptable, or should the offset use I*(1 − x)?" },
    { id: "q-hbr", ch: "stage1", q: "Are the interim abundance maps on an absolute human-biting-rate scale, and against which catch method were they calibrated?" },
    { id: "q-ir", ch: "interventions", q: "Is the Symons 0.46 retained-effect estimate acceptable as an interim mapping from IR-cube susceptibility to net effect, before a hut-trial mapping?" },
    { id: "q-survival", ch: "stage1", q: "Can I have the dehydrated g(T, H) survival function from the An. stephensi work, and the humidity biting lookup when it is ready?" },
    { id: "q-window", ch: "likelihood", q: "With six-month infections, should the detectability kernel q follow the infection's duration (q(d) ∝ e^{−rd}) instead of a 30-day window, so survey prevalence matches the ODE?" },
    { id: "q-params", ch: "stage1", q: "The simulation now uses literature values (r = 1/180, b = c = 0.5, EIP 10 days, m = 0.2), and the case likelihood multiplies by days per month. Agree before rerunning the study?" },
    { id: "q-centre", ch: "stage2", q: "Should α-centring be the default everywhere? It changes α to α + mean(ε)." }
  ],

  formulas: [
    { group: "Stage 1", items: [
      ["Human infection", "dx/dt = m a b z (1 − x) − r x"],
      ["Mosquito infection", "dz/dt = a c x (1 − z) − g z"],
      ["Mechanistic incidence rate", "I* = m · a · b · z", "compute_mechanistic_prediction()"],
      ["Basic reproduction number", "R₀ = m a² b c S / (g r),  S = survival through the EIP"],
      ["Equilibrium prevalence", "x* = (R₀ − 1) / (R₀ + a c / g)"],
      ["ITN effect on m", "m (1 − n′ k),  n′ = n (1 − 0.46 (1 − q))", "apply_interventions()"],
      ["Surviving the EIP", "(1 + g n / 4)<sup>−4</sup> ≈ e<sup>−g n</sup>", "eip_survival()"]
    ]},
    { group: "Stage 2", items: [
      ["Latent incidence", "log I = α + log I* + ε", "fit_epiwave_gp()"],
      ["Spatial kernel", "k(d) = σ² (1 + √5 d/φ + 5d²/3φ²) e<sup>−√5 d/φ</sup>", "mat52()"],
      ["AR(1) in time", "ε<sub>t</sub> = θ ε<sub>t−1</sub> + f<sub>t</sub>", "ar1()"],
      ["Stationary variance", "τ² = σ² / (1 − θ²)", "PRIORS$tau2"]
    ]},
    { group: "Likelihoods", items: [
      ["Cases", "C ~ Poisson(γ · I · N · 30)  (I per day, cases per month)"],
      ["Detectable infections", "Λ<sub>t</sub> = Σ<sub>k</sub> I<sub>t−k</sub> w<sub>k</sub>", "build_detectability_matrix()"],
      ["Survey prevalence", "p = 1 − e<sup>−Λ</sup>"],
      ["Survey positives", "Y ~ Binomial(T, p)"]
    ]},
    { group: "Priors", items: [
      ["Intercept", "α ~ Normal(0, 1)"],
      ["Reporting rate", "γ ~ Normal(0.1, 0.05), truncated at 0.001"],
      ["Stationary variance", "τ² ~ LogNormal(−0.5, 0.5)"],
      ["AR(1) correlation", "θ ~ Beta(2, 2)"],
      ["Lengthscale", "φ ~ LogNormal(0.5, 0.5)"]
    ]}
  ]
};
