module sphharm

using AssociatedLegendrePolynomials, LinearAlgebra, Tullio
using ProgressMeter, OffsetArrays


# naive implementation, it takes way more time then GLQ-based approach in SHTOOLS 
function direct_analysis(grid, values, N)
	t = eltype(values)
    Cnk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    Snk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    areas = prod(grid.step) * cos.(grid.second_axis)
    @showprogress for n in 0:N
        for k in 0:n
            s₁ = zero(t)
            s₂ = zero(t)
            for i in eachindex(grid.second_axis)
				pnk = correct_lgn(sin(grid.second_axis[i]), n, k)
				area = areas[i]
				for j in eachindex(grid.first_axis)
                	sin_kλ, cos_kλ = sincos(k * grid.first_axis[j])
                	s₁ += values[i,j] * pnk * cos_kλ * area
                	s₂ += values[i,j] * pnk * sin_kλ * area
            	end
			end
        Cnk[n,k] = s₁
        Snk[n,k] = s₂
        end
    end
    return 1/4π * Cnk.parent, 1/4π * Snk.parent
end



function optimized_analysis(grid, values, N)
	t = eltype(values)
    weighted_values = values * prod(grid.step) .* cos.(grid.second_axis)
    cos_kλ = compute_order_series(grid.first_axis, N, f = cos)
    sin_kλ = compute_order_series(grid.first_axis, N, f = sin)
    
    # preallocate output
    pnk = zeros(N+1, N+1)
    # Cnk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    # Snk = OffsetArray(zeros(t, N+1, N+1), 0:N, 0:N)
    Cnk = zeros(t, N+1, N+1)
    Snk = zeros(t, N+1, N+1)

    @showprogress for (i,φ) in enumerate(grid.second_axis)
        correct_lgn!(pnk, sin(φ), N)
        row_values = view(weighted_values, i, :)
        for n in 0:N
            for k in 0:n 
                Cnk[n+1,k+1] += pnk[n+1,k+1] * dot(row_values, view(cos_kλ, k+1, :))
                Snk[n+1,k+1] += pnk[n+1,k+1] * dot(row_values, view(sin_kλ, k+1, :))
            end
        end
    end
    return 1/4π * Cnk, 1/4π * Snk
end



# must fix the normalization
function correct_lgn(x, n, k)
    pnk = legendre(LegendreFourPiNorm(), n, k, x)
    if k == 0
        return pnk
    else
        if isodd(k)
            return -pnk * √2 
        else
            return pnk * √2
        end
    end
end


# the same
function correct_lgn!(out, x, N)
	legendre!(LegendreFourPiNorm(), out, N, N, x)
    @views for i in 2:N+1
		v = out[i:end,i]
        if isodd(i)
            map!(x -> √2 * x, out[i:end,i], v)
        else
            map!(x -> -√2 * x, out[i:end,i], v)
        end
    return out
    end
end


function lumped_coefficients!(aₖ, bₖ, pnk, cnk, snk, N)
    for k in 0:N
		idx = k + 1
        @views pnk_vec = pnk[idx:end,idx] 
        @views aₖ[idx] = dot(pnk_vec, cnk[idx:end,idx])
        @views bₖ[idx] = dot(pnk_vec, snk[idx:end,idx])
    end
    return aₖ, bₖ
end


function compute_order_series(vector, order_max; f)
    return f.(vector' .* collect(0:order_max))
end


function sum_tullio!(out, aₖ, cos_kλ, bₖ, sin_kλ)
	@tullio threads = false out[j] = aₖ[n] * cos_kλ[n,j] + bₖ[n] * sin_kλ[n,j]
end	


function synthesis!(out, grid, cnk, snk, N)
    cos_kλ = compute_order_series(grid.first_axis, N, f = cos)
    sin_kλ = compute_order_series(grid.first_axis, N, f = sin)
    pnk = zeros(N+1,N+1)
	aₖ = zeros(N+1)
    bₖ = zeros(N+1)
    @showprogress for (i,φ) in enumerate(grid.second_axis)
        correct_lgn!(pnk, sin(φ), N)
        lumped_coefficients!(aₖ, bₖ, pnk, cnk, snk, N)
        sum_tullio!(view(out, i, :), aₖ, cos_kλ, bₖ, sin_kλ)
    end
    return out
end


function synthesis(point, delta_factor, cnk, snk, N)
    λ,φ = point
    f = -6 / (5.5134e+3 * 6371000)
    cos_kλ = compute_order_series(λ, N, f = cos)
    sin_kλ = compute_order_series(λ, N, f = sin)

    pnk = zeros(N+1,N+1)
    correct_lgn!(pnk, sin(φ), N)
    s = sum((cnk .* cos_kλ' + snk .* sin_kλ') .* pnk, dims = 2)
    
    return f * dot(delta_factor, s)

end


end