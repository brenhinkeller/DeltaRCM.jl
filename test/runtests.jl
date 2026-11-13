using DeltaRCM
using Test
using NCDatasets
using JLD2
using Random

const SMALL = (; Length = 1000, Width = 2000, dx = 50, Np_water = 200, Np_sed = 200, seed = 11)

@testset "DeltaRCM.jl" begin
    @testset "parameters" begin
        p = DeltaParameters()
        @test p.Length == 5000 && p.Np_water == 2000 && p.seed === nothing
        @test DeltaParameters(Np_water = 10).Np_water == 10
        @test_throws ArgumentError DeltaParameters(Dict(:Np_water => 1.5))
        @test_throws ArgumentError DeltaParameters(Dict(:f_bedload => 2))
        @test_logs (:warn, r"Unused") DeltaParameters(Dict(:not_a_param => 1))
        @test DeltaParameters(Dict("hb" => nothing, "seed" => 3)).seed == 3
        mktempdir() do dir
            path = joinpath(dir, "c.yaml")
            write(path, "Length: 1000.\nWidth: 2000.\nseed: 5\nS0: 2e-4\nhb: null\n")
            q = DeltaParameters(path)
            @test q.Length == 1000 && q.seed == 5 && q.S0 == 2e-4 && q.hb === nothing
            @test_logs (:warn, r"also specified") DeltaParameters(path; seed = 6)
            @test DeltaParameters(path; Np_water = 7).Np_water == 7
        end
        @test DeltaRCM.with(p; dx = 25.0).dx == 25
    end

    @testset "derived constants and rounding" begin
        m = @test_logs (:warn, r"Length") (:warn, r"L0_meters") DeltaModel(; Length = 1010, Width = 2000, dx = 50,
                                                                           L0_meters = 160)
        @test size(m) == (20, 40)
        @test m.L0 == 3 && m.N0 == 5
        @test m.dry_depth ≈ 0.1
        @test m.dt ≈ m.dVs / m.Qs0
        @test m.gamma ≈ 9.81 * 0.00015 * 50 / 1
        # default domain
        d = DeltaModel()
        @test size(d) == (100, 200)
        @test d.CTR == 100 && d.L0 == 3 && d.N0 == 5
    end

    @testset "initial domain" begin
        m = DeltaModel(; SMALL...)
        ct = m.cell_type
        @test all(ct[1:m.L0, :] .∈ Ref((DeltaRCM.WALL, DeltaRCM.CHANNEL)))
        @test count(==(DeltaRCM.CHANNEL), ct) == m.L0 * m.N0
        @test length(m.inlet) == m.N0
        @test all(ct[end, :] .== DeltaRCM.EDGE)
        @test all(ct[m.L0+1:end, 1] .== DeltaRCM.EDGE) && all(ct[m.L0+1:end, end] .== DeltaRCM.EDGE)
        @test all(m.depth[ct.==DeltaRCM.CHANNEL] .== m.h0)
        @test m.eta == m.stage .- m.depth
        @test all(DeltaRCM.interior(m.f.cell_type) .== ct)
        # halo is wall
        @test all(m.f.cell_type[1, :] .== DeltaRCM.WALL) && all(m.f.cell_type[:, end] .== DeltaRCM.WALL)
    end

    @testset "water weights" begin
        m = DeltaModel(; SMALL...)
        f, c = m.f, m.c
        # a cell in the channel mouth: row L0 (the last channel row), centre column
        idx = (c.L0 + 1) + (c.CTR + 1 - 1) * c.Lp
        w, tot = DeltaRCM.water_weights(f, c, idx)
        @test w[5] == 0 && tot ≈ sum(w) && tot > 0
        # flow is downstream (+i) so the downstream neighbour (k = 6) outweighs the upstream one (k = 4)
        @test w[6] > w[4]
        # cells never step into walls: at the inlet row the upstream neighbours are walls/halo
        top = 2 + (c.CTR + 1 - 1) * c.Lp
        wt, tt = DeltaRCM.water_weights(f, c, top)
        @test wt[1] == 0 && wt[4] == 0 && wt[7] == 0 && tt > 0
    end

    @testset "reproducibility" begin
        a = DeltaModel(; SMALL...)
        b = DeltaModel(; SMALL...)
        run!(a, 5); run!(b, 5)
        @test a.eta == b.eta && a.stage == b.stage && a.qx == b.qx
        c = DeltaModel(; SMALL..., seed = 12)
        run!(c, 5)
        @test c.eta != a.eta
    end

    @testset "mass conservation" begin
        m = DeltaModel(; SMALL...)
        for _ in 1:15
            e0 = copy(m.eta)
            step!(m)
            Vin = (floor(Int, m.Np_sed * m.f_bedload) + floor(Int, m.Np_sed * (1 - m.f_bedload))) * m.Vp_sed
            t = m.tally
            out = t.exported_sand + t.exported_mud + t.lost_sand + t.lost_mud
            dV = sum(m.eta .- e0) * m.dx^2
            @test dV ≈ Vin - out + m.Vp_inletbc rtol = 1e-9
        end
    end

    @testset "physical sanity" begin
        m = DeltaModel(; SMALL...)
        run!(m, 20)
        @test all(isfinite, m.eta) && all(isfinite, m.stage) && all(isfinite, m.qx) && all(isfinite, m.uw)
        @test all(m.depth .>= 0)
        @test all(m.uw .<= m.u_max + sqrt(eps()))
        @test all(m.stage .>= m.H_SL .- sqrt(eps()))
        @test maximum(m.eta) > m.eta_init[end, end] + 0.1   # delta has built up somewhere
        @test m.time ≈ 20 * m.dt && m.time_iter == 20
        @test all((0 .<= m.sand_frac) .& (m.sand_frac .<= 1))
    end

    @testset "Float32" begin
        m = DeltaModel{Float32}(; SMALL...)
        run!(m, 5)
        @test eltype(m.eta) == Float32 && all(isfinite, m.eta)
    end

    @testset "sea level rise and subsidence" begin
        m = DeltaModel(; SMALL..., SLR = 1e-7, toggle_subsidence = true, subsidence_rate = 1e-7,
                       start_subsidence = 0.0)
        run!(m, 3)
        @test m.H_SL ≈ 3 * 1e-7 * m.dt
        @test m.sigma[end, end] ≈ 1e-7 * m.dt && all(m.sigma[1:m.L0, :] .== 0)
    end

    @testset "finalize" begin
        m = DeltaModel(; SMALL...)
        finalize!(m)
        @test_throws ErrorException step!(m)
        @test_throws ErrorException finalize!(m)
    end

    @testset "NetCDF output (pyDeltaRCM layout)" begin
        mktempdir() do dir
            m = DeltaModel(; SMALL..., out_dir = dir, save_eta_grids = true, save_depth_grids = true,
                           save_velocity_components = true, save_dt = 1.0)
            run!(m, 3)
            finalize!(m)
            NCDataset(joinpath(dir, "pyDeltaRCM_output.nc")) do ds
                @test ds.attrib["sandsuet_version"] == "1.0.0"
                @test Set(keys(ds)) >= Set(["eta", "depth", "velocity_x", "velocity_y", "seconds", "x", "y"])
                @test dimnames(ds["eta"]) == ("y", "x", "seconds")       # (seconds, x, y) in NetCDF/Python order
                @test size(ds["eta"]) == (40, 20, 4)
                @test ds["seconds"][:] ≈ (0:3) .* m.dt
                @test ds["eta"][:, :, 4] ≈ permutedims(m.eta)
                @test ds["eta"][:, :, 1] ≈ permutedims(m.eta_init)
                aux = ds.group["auxdata"]
                @test aux["CTR"][] == m.CTR - 1 && aux["L0"][] == 3 && aux["dx"][] == 50
                @test aux["cell_type"][:, :] == permutedims(m.cell_type)
                @test length(aux["H_SL"]) == 4
                @test ds["eta"].attrib["long_name"] == "channel_bottom__elevation"
            end
            # existing file is protected unless clobbering
            @test_throws ErrorException DeltaModel(; SMALL..., out_dir = dir, save_eta_grids = true)
            m2 = @test_logs (:warn, r"Replacing") DeltaModel(; SMALL..., out_dir = dir, save_eta_grids = true,
                                                           clobber_netcdf = true)
            finalize!(m2)
            # save_dt thins the output
            m3 = DeltaModel(; SMALL..., out_dir = joinpath(dir, "b"), save_eta_grids = true, save_dt = 2.5 * 25000)
            run!(m3, 5); finalize!(m3)
            NCDataset(joinpath(dir, "b", "pyDeltaRCM_output.nc")) do ds
                # saves when >= save_dt has elapsed: t = 0, then after step 3 (75000 s >= 62500 s)
                @test ds["seconds"][:] ≈ [0, 3 * m3.dt]
            end
        end
    end

    @testset "checkpoint and resume" begin
        mktempdir() do dir
            cfg = (; SMALL..., out_dir = dir, save_checkpoint = true, save_eta_grids = true, save_dt = 1.0)
            ref = DeltaModel(; SMALL...)
            run!(ref, 6)
            a = DeltaModel(; cfg...)
            run!(a, 3)
            finalize!(a)
            @test isfile(joinpath(dir, "checkpoint.jld2"))
            b = DeltaModel(; cfg..., resume_checkpoint = true)
            @test b.iter == 3 && b.time ≈ 3 * b.dt
            run!(b, 3)
            @test b.eta == ref.eta && b.stage == ref.stage && b.sand_frac == ref.sand_frac
            finalize!(b)
            NCDataset(joinpath(dir, "pyDeltaRCM_output.nc")) do ds
                @test length(ds["seconds"]) == 7             # t = 0 plus six steps, no duplicates
                @test ds["eta"][:, :, 7] ≈ permutedims(ref.eta)
            end
        end
        @test_throws ArgumentError DeltaModel(; SMALL..., out_dir = mktempdir(), resume_checkpoint = true)
    end

    @testset "run durations" begin
        m = DeltaModel(; SMALL...)
        run!(m; timesteps = 2)
        @test m.iter == 2
        run!(m; time = 1.5 * m.dt)                   # one whole step plus a half step
        @test m.iter == 4 && m.time ≈ 3.5 * m.dt
        run!(m; time_years = 0.0)
        @test m.iter == 4
        @test_throws ArgumentError run!(m)
        @test_throws ArgumentError run!(m; time = 1.0, timesteps = 1)
        # a shortened step carries a proportional sediment load: mass is still conserved
        m2 = DeltaModel(; SMALL...)
        e0 = copy(m2.eta)
        step!(m2; dt = 0.5 * m2.dt)
        t = m2.tally
        dV = sum(m2.eta .- e0) * m2.dx^2
        Vin = (floor(Int, m2.Np_sed * m2.f_bedload) + floor(Int, m2.Np_sed * (1 - m2.f_bedload))) * m2.Vp_sed / 2
        @test dV ≈ Vin - (t.exported_sand + t.exported_mud + t.lost_sand + t.lost_mud) + m2.Vp_inletbc rtol = 1e-9
        @test m2.time ≈ 0.5 * m2.dt && m2.c.dt == m2.dt0       # nominal constants restored
    end

    @testset "run_model and jobs" begin
        cfg = Dict(:Length => 1000, :Width => 2000, :dx => 50, :Np_water => 100, :Np_sed => 100, :seed => 4)
        m = run_model(cfg; timesteps = 2)
        @test m.iter == 2 && m.is_finalized
        @test_throws ArgumentError run_model(cfg)                           # no duration
        mktempdir() do dir
            yml = joinpath(dir, "c.yaml")
            write(yml, "Length: 1000\nWidth: 2000\ndx: 50\nNp_water: 100\nNp_sed: 100\nseed: 4\ntimesteps: 2\n")
            @test run_model(yml).eta == m.eta
        end
        # expansion
        cfg2 = merge(cfg, Dict(:out_dir => "root", :timesteps => 1))
        delete!(cfg2, :seed)
        jobs = expand_jobs(merge(cfg2, Dict(:matrix => Dict(:f_bedload => [0.3, 0.5], :u0 => [0.8, 1.0, 1.2]))))
        @test length(jobs) == 6 && jobs[1][:out_dir] == joinpath("root", "job_000")
        @test [j[:u0] for j in jobs[1:3]] == [0.8, 1.0, 1.2] && jobs[4][:f_bedload] == 0.5
        @test length(expand_jobs(merge(cfg2, Dict(:ensemble => 4)))) == 4
        sets = expand_jobs(merge(cfg2, Dict(:set => [Dict(:u0 => 1.0, :h0 => 4.0), Dict(:u0 => 1.2, :h0 => 5.0)], :ensemble => 2)))
        @test length(sets) == 4 && length(unique(j[:seed] for j in sets)) == 4
        @test_throws ArgumentError expand_jobs(merge(cfg2, Dict(:set => [], :matrix => Dict())))
        @test_throws ArgumentError expand_jobs(Dict(:matrix => Dict(:u0 => [1, 2])))      # needs out_dir
        mktempdir() do dir
            cfg3 = merge(cfg, Dict(:out_dir => joinpath(dir, "sweep"), :timesteps => 2, :ensemble => 3))
            delete!(cfg3, :seed)
            vols = run_jobs(m -> sum(m.eta), cfg3; rng = Random.Xoshiro(1))
            @test length(vols) == 3 && length(unique(vols)) == 3
            @test readlines(joinpath(dir, "sweep", "jobs_parameters.txt"))[1] == "job_id, seed"
            @test_throws ArgumentError run_jobs(cfg3)                    # directory exists
            par = run_jobs(m -> sum(m.eta), merge(cfg3, Dict(:out_dir => joinpath(dir, "par"))),
                           rng = Random.Xoshiro(1), parallel = true)
            @test par == vols                                                         # same seeds, same results
        end
    end

    @testset "hooks" begin
        log = Symbol[]
        hooks = NamedTuple{DeltaRCM.HOOK_POINTS}(ntuple(i -> (m -> push!(log, DeltaRCM.HOOK_POINTS[i])), 7))
        m = DeltaModel(; SMALL..., hooks)
        step!(m)
        @test log == collect(DeltaRCM.HOOK_POINTS)
        @test_throws ArgumentError DeltaModel(; SMALL..., hooks = (; nope = identity))
        # a hook that makes the left half of the basin unattractive to water changes the outcome
        base = DeltaModel(; SMALL...); run!(base, 5)
        w = DeltaModel(; SMALL..., hooks = (; init = m -> (m.mod_water_weight[:, 1:m.CTR] .= 0.2)))
        run!(w, 5)
        @test w.eta != base.eta
    end

    @testset "weight tables" begin
        m = DeltaModel(; SMALL...)
        run!(m, 6)
        f, c, sc = m.f, m.c, m.scratch
        # water: cached weights equal a fresh computation, after parcels have run an iteration
        DeltaRCM.run_water_parcels!(sc, f, c, m.inlet, m.seed, m.iter, 1, m.H_SL)
        visited = findall(vec(sc.wstamp) .== sc.wgen[])
        @test length(visited) > 50
        for idx in visited[1:10:end]
            w, t = DeltaRCM.water_weights(f, c, idx)
            @test sc.wtab[idx] == w && sc.wtot[idx] == t
        end
        # sediment: the attractiveness table stays consistent with the bed after each phase
        for kind in (DeltaRCM.Sand(), DeltaRCM.Mud())
            DeltaRCM.route_sediment_parcels!(kind, f, c, m.inlet, m.seed, m.iter, m.tally)
            fresh = [DeltaRCM.sediment_attract(kind, f, c, n) for n in eachindex(f.depth)]
            @test vec(f.sed_attract) == fresh
        end
    end
end
