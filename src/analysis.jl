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


