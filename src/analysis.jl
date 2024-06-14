function direct_analysis(grid, values, N)
    Cnk = OffsetArray(zeros(N+1, N+1), 0:N, 0:N)
    Snk = OffsetArray(zeros(N+1, N+1), 0:N, 0:N)
    areas = @. grid.step[1] * grid.step[2] * cos(grid.second_axis)
    # axis_vec = view(grid.first_axis, 1, :)
    #cos_kλ = compute_order_series(axis_vec, N, f = cos)
    #sin_kλ = compute_order_series(axis_vec, N, f = sin)
    for n in tqdm(0:N)
        for k in 0:n
            s₁ = 0.0
            s₂ = 0.0
            for i in eachindex(values)
                λ,φ = grid.first_axis[i], grid.second_axis[i]
                pnk = correct_lgn(sin(φ), n, k)
                sin_kλ, cos_kλ = sincos(k * λ)
                s₁ += values[i] * pnk * cos_kλ * areas[i]
                s₂ += values[i] * pnk * sin_kλ * areas[i]
            end
            Cnk[n,k] = 1/4π * s₁
            Snk[n,k] = 1/4π * s₂
        end
    end
    return Cnk, Snk
end


