module Takums

import Base: AbstractFloat, Int, Int8, Int16, Int32, Int64, Integer, Signed,
	Unsigned, reinterpret

import Printf

using libtakum_jll

export Takum8, Takum16, Takum32, Takum64, AnyTakum,
       NaR8, NaR16, NaR32, NaR64, takum, isnar

# type definitions
primitive type Takum8  <: AbstractFloat 8 end
primitive type Takum16 <: AbstractFloat 16 end
primitive type Takum32 <: AbstractFloat 32 end
primitive type Takum64 <: AbstractFloat 64 end
AnyTakum = Union{Takum8, Takum16, Takum32, Takum64}

# NaR representations
const NaR8  = Base.bitcast(Takum8,  0x80)
const NaR16 = Base.bitcast(Takum16, 0x8000)
const NaR32 = Base.bitcast(Takum32, 0x80000000)
const NaR64 = Base.bitcast(Takum64, 0x8000000000000000)

# array takum types with their corresponding integer types for
# metaprogramming within the module
takum_types = [
	(:Takum8, :Int8),
	(:Takum16, :Int16),
	(:Takum32, :Int32),
	(:Takum64, :Int64),
]

# integer reinterpret casts
Base.reinterpret(::Type{Unsigned}, x::Takum8)  = reinterpret(UInt8, x)
Base.reinterpret(::Type{Unsigned}, x::Takum16) = reinterpret(UInt16, x)
Base.reinterpret(::Type{Unsigned}, x::Takum32) = reinterpret(UInt32, x)
Base.reinterpret(::Type{Unsigned}, x::Takum64) = reinterpret(UInt64, x)

Base.reinterpret(::Type{Signed}, x::Takum8)  = reinterpret(Int8, x)
Base.reinterpret(::Type{Signed}, x::Takum16) = reinterpret(Int16, x)
Base.reinterpret(::Type{Signed}, x::Takum32) = reinterpret(Int32, x)
Base.reinterpret(::Type{Signed}, x::Takum64) = reinterpret(Int64, x)

Base.uinttype(::Type{Takum8})  = UInt8
Base.uinttype(::Type{Takum16}) = UInt16
Base.uinttype(::Type{Takum32}) = UInt32
Base.uinttype(::Type{Takum64}) = UInt64

# the only floating-point property that makes sense to implement is signbit()
Base.signbit(t::Takum8)  = (reinterpret(Unsigned, t) & 0x80) !== 0x00
Base.signbit(t::Takum16) = (reinterpret(Unsigned, t) & 0x8000) !== 0x0000
Base.signbit(t::Takum32) = (reinterpret(Unsigned, t) & 0x80000000) !== 0x00000000
Base.signbit(t::Takum64) = (reinterpret(Unsigned, t) & 0x8000000000000000) !== 0x0000000000000000

# left undefined are sign_mask, exponent_mask, exponent_one, exponent_half,
# significand_mask, exponent_bias, exponent_bits, significand_bits, significand,
# exponent, decompose, frexp, ldexp

Base.iszero(t::Takum8)  = reinterpret(Unsigned, t) === 0x00
Base.iszero(t::Takum16) = reinterpret(Unsigned, t) === 0x0000
Base.iszero(t::Takum32) = reinterpret(Unsigned, t) === 0x00000000
Base.iszero(t::Takum64) = reinterpret(Unsigned, t) === 0x0000000000000000

Base.isone(t::Takum8)  = reinterpret(Unsigned, t) === 0x40
Base.isone(t::Takum16) = reinterpret(Unsigned, t) === 0x4000
Base.isone(t::Takum32) = reinterpret(Unsigned, t) === 0x40000000
Base.isone(t::Takum64) = reinterpret(Unsigned, t) === 0x4000000000000000

Base.isfinite(t::Takum8)  = reinterpret(Unsigned, t) !== 0x80
Base.isfinite(t::Takum16) = reinterpret(Unsigned, t) !== 0x8000
Base.isfinite(t::Takum32) = reinterpret(Unsigned, t) !== 0x80000000
Base.isfinite(t::Takum64) = reinterpret(Unsigned, t) !== 0x8000000000000000

Base.isnan(t::Takum8)  = reinterpret(Unsigned, t) === 0x80
Base.isnan(t::Takum16) = reinterpret(Unsigned, t) === 0x8000
Base.isnan(t::Takum32) = reinterpret(Unsigned, t) === 0x80000000
Base.isnan(t::Takum64) = reinterpret(Unsigned, t) === 0x8000000000000000

isnar(t::Takum8)  = Base.isnan(t)
isnar(t::Takum16) = Base.isnan(t)
isnar(t::Takum32) = Base.isnan(t)
isnar(t::Takum64) = Base.isnan(t)

Base.issubnormal(t::AnyTakum) = false
Base.ispow2(t::AnyTakum) = Base.isone(t)
Base.iseven(t::AnyTakum) = Base.iszero(t)
Base.isodd(t::AnyTakum) = Base.isone(t) || Base.isone(-t)

# precision
_mantissa_bit_count(t::Takum8)  = @ccall libtakum.takum8_precision(reinterpret(Signed, t)::Int8)::UInt8
_mantissa_bit_count(t::Takum16) = @ccall libtakum.takum16_precision(reinterpret(Signed, t)::Int16)::UInt8
_mantissa_bit_count(t::Takum32) = @ccall libtakum.takum32_precision(reinterpret(Signed, t)::Int32)::UInt8
_mantissa_bit_count(t::Takum64) = @ccall libtakum.takum64_precision(reinterpret(Signed, t)::Int64)::UInt8

function Base.precision(t::AnyTakum; base::Integer = 2)
	base > 1 || throw(DomainError(base, "`base` cannot be less than 2."))
	m = _mantissa_bit_count(t)
	return base == 2 ? Int(m) : floor(Int, m / log2(base))
end

# For the types we determine the precision of the zero of said type, as this returns
# the worst-case precision, consistent with what you obtain with the respective
# IEEE 754 precision functions
Base.precision(T::Type{<:AnyTakum}; base::Integer = 2) = precision(zero(T); base)

# eps (follow definition; for the types it is simply eps(Takum_(1.0))
Base.eps(t::AnyTakum)     = max(t - prevfloat(t), nextfloat(t) - t)
Base.eps(::Type{Takum8})  = Base.bitcast(Takum8,  0x2b)
Base.eps(::Type{Takum16}) = Base.bitcast(Takum16, 0x1f2f)
Base.eps(::Type{Takum32}) = Base.bitcast(Takum32, 0x160bc2b0)
Base.eps(::Type{Takum64}) = Base.bitcast(Takum64, 0x0d7a50987ab4d7df)

# rounding
Base.round(t::AnyTakum) = typeof(t)(Base.round(Float64(t)))

Base.round(t::AnyTakum, r::RoundingMode{:ToZero})  = typeof(t)(Base.trunc(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Down})    = typeof(t)(Base.floor(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Up})      = typeof(t)(Base.ceil(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Nearest}) = typeof(t)(Base.round(Float64(t)))

Base.trunc(t::AnyTakum) = Base.signbit(t) ? Base.ceil(t) : Base.floor(t)
Base.trunc(::Type{T}, t::AnyTakum) where {T <: Integer} = Base.trunc(T, Float64(t))

# type limits
Base.typemin(::Type{Takum8})  = NaR8
Base.typemin(::Type{Takum16}) = NaR16
Base.typemin(::Type{Takum32}) = NaR32
Base.typemin(::Type{Takum64}) = NaR64

Base.typemax(::Type{Takum8})  = Base.bitcast(Takum8,  0x7f)
Base.typemax(::Type{Takum16}) = Base.bitcast(Takum16, 0x7fff)
Base.typemax(::Type{Takum32}) = Base.bitcast(Takum32, 0x7fffffff)
Base.typemax(::Type{Takum64}) = Base.bitcast(Takum64, 0x7fffffffffffffff)

Base.floatmin(::Type{Takum8})  = Base.bitcast(Takum8,  0x01)
Base.floatmin(::Type{Takum16}) = Base.bitcast(Takum16, 0x0001)
Base.floatmin(::Type{Takum32}) = Base.bitcast(Takum32, 0x00000001)
Base.floatmin(::Type{Takum64}) = Base.bitcast(Takum64, 0x0000000000000001)

Base.floatmax(::Type{Takum8})  = Base.bitcast(Takum8,  0x7f)
Base.floatmax(::Type{Takum16}) = Base.bitcast(Takum16, 0x7fff)
Base.floatmax(::Type{Takum32}) = Base.bitcast(Takum32, 0x7fffffff)
Base.floatmax(::Type{Takum64}) = Base.bitcast(Takum64, 0x7fffffffffffffff)

Base.maxintfloat(::Type{T}) where {T<:AnyTakum} = T(1.0)

# conversions from floating-point
Takum8(f::Float16) = Base.bitcast(Takum8, @ccall libtakum.takum8_from_float32(Float32(f)::Float32)::Int8)
Takum8(f::Float32) = Base.bitcast(Takum8, @ccall libtakum.takum8_from_float32(f::Float32)::Int8)
Takum8(f::Float64) = Base.bitcast(Takum8, @ccall libtakum.takum8_from_float64(f::Float64)::Int8)

Takum16(f::Float16) = Base.bitcast(Takum16, @ccall libtakum.takum16_from_float32(Float32(f)::Float32)::Int16)
Takum16(f::Float32) = Base.bitcast(Takum16, @ccall libtakum.takum16_from_float32(f::Float32)::Int16)
Takum16(f::Float64) = Base.bitcast(Takum16, @ccall libtakum.takum16_from_float64(f::Float64)::Int16)

Takum32(f::Float16) = Base.bitcast(Takum32, @ccall libtakum.takum32_from_float32(Float32(f)::Float32)::Int32)
Takum32(f::Float32) = Base.bitcast(Takum32, @ccall libtakum.takum32_from_float32(f::Float32)::Int32)
Takum32(f::Float64) = Base.bitcast(Takum32, @ccall libtakum.takum32_from_float64(f::Float64)::Int32)

Takum64(f::Float16) = Base.bitcast(Takum64, @ccall libtakum.takum64_from_float32(Float32(f)::Float32)::Int64)
Takum64(f::Float32) = Base.bitcast(Takum64, @ccall libtakum.takum64_from_float32(f::Float32)::Int64)
Takum64(f::Float64) = Base.bitcast(Takum64, @ccall libtakum.takum64_from_float64(f::Float64)::Int64)

# conversion from integers with promote rules
Takum8(i::Integer)  = Base.convert(Takum8,  Base.convert(Float64, i))
Takum16(i::Integer) = Base.convert(Takum16, Base.convert(Float64, i))
Takum32(i::Integer) = Base.convert(Takum32, Base.convert(Float64, i))
Takum64(i::Integer) = Base.convert(Takum64, Base.convert(Float64, i))
Base.promote_rule(T::Type{<:AnyTakum}, ::Type{<:Integer}) = T

# conversions to floating-point
Base.Float16(t::Takum8)  = Float16(@ccall libtakum.takum8_to_float32(reinterpret(Signed, t)::Int8)::Float32)
Base.Float16(t::Takum16) = Float16(@ccall libtakum.takum16_to_float32(reinterpret(Signed, t)::Int16)::Float32)
Base.Float16(t::Takum32) = Float16(@ccall libtakum.takum32_to_float32(reinterpret(Signed, t)::Int32)::Float32)
Base.Float16(t::Takum64) = Float16(@ccall libtakum.takum64_to_float32(reinterpret(Signed, t)::Int64)::Float32)
Base.Float32(t::Takum8)  = @ccall libtakum.takum8_to_float32(reinterpret(Signed, t)::Int8)::Float32
Base.Float32(t::Takum16) = @ccall libtakum.takum16_to_float32(reinterpret(Signed, t)::Int16)::Float32
Base.Float32(t::Takum32) = @ccall libtakum.takum32_to_float32(reinterpret(Signed, t)::Int32)::Float32
Base.Float32(t::Takum64) = @ccall libtakum.takum64_to_float32(reinterpret(Signed, t)::Int64)::Float32
Base.Float64(t::Takum8)  = @ccall libtakum.takum8_to_float64(reinterpret(Signed, t)::Int8)::Float64
Base.Float64(t::Takum16) = @ccall libtakum.takum16_to_float64(reinterpret(Signed, t)::Int16)::Float64
Base.Float64(t::Takum32) = @ccall libtakum.takum32_to_float64(reinterpret(Signed, t)::Int32)::Float64
Base.Float64(t::Takum64) = @ccall libtakum.takum64_to_float64(reinterpret(Signed, t)::Int64)::Float64

# conversion to integer
Base.unsafe_trunc(T::Type{<:Integer}, t::AnyTakum)  = Base.unsafe_trunc(T, Float64(t))
Base.round(I::Type{<:Integer}, t::AnyTakum) = Base.round(I, Float64(t))

function (::Type{I})(t::AnyTakum) where I <: Integer
	if t == -1
		return I(-1)
	elseif t == 0
		return I(0)
	elseif t == 1
		return I(1)
	else
		throw(InexactError(:round, I, t))
	end
end

# inter-takum conversions
Takum8(t::Takum16) = Base.bitcast(Takum8, UInt8(reinterpret(Unsigned, t) >> 8) + ((reinterpret(Unsigned, t) & (UInt16(1) << 7)) != 0))
Takum8(t::Takum32) = Base.bitcast(Takum8, UInt8(reinterpret(Unsigned, t) >> 24) + ((reinterpret(Unsigned, t) & (UInt32(1) << 23)) != 0))
Takum8(t::Takum64) = Base.bitcast(Takum8, UInt8(reinterpret(Unsigned, t) >> 56) + ((reinterpret(Unsigned, t) & (UInt64(1) << 55)) != 0))

Takum16(t::Takum8)  = Base.bitcast(Takum16, UInt16(reinterpret(Unsigned, t)) << 8)
Takum16(t::Takum32) = Base.bitcast(Takum16, UInt16(reinterpret(Unsigned, t) >> 16) + ((reinterpret(Unsigned, t) & (UInt32(1) << 15)) != 0))
Takum16(t::Takum64) = Base.bitcast(Takum16, UInt16(reinterpret(Unsigned, t) >> 48) + ((reinterpret(Unsigned, t) & (UInt32(1) << 47)) != 0))

Takum32(t::Takum8)  = Base.bitcast(Takum32, UInt32(reinterpret(Unsigned, t)) << 24)
Takum32(t::Takum16) = Base.bitcast(Takum32, UInt32(reinterpret(Unsigned, t)) << 16)
Takum32(t::Takum64) = Base.bitcast(Takum32, UInt32(reinterpret(Unsigned, t) >> 32) + ((reinterpret(Unsigned, t) & (UInt64(1) << 31)) != 0))

Takum64(t::Takum8)  = Base.bitcast(Takum64, UInt64(reinterpret(Unsigned, t)) << 56)
Takum64(t::Takum16) = Base.bitcast(Takum64, UInt64(reinterpret(Unsigned, t)) << 48)
Takum64(t::Takum32) = Base.bitcast(Takum64, UInt64(reinterpret(Unsigned, t)) << 32)

# inter-takum promote rules
Base.promote_rule(::Type{Takum16}, ::Type{Takum8}) = Takum16

Base.promote_rule(::Type{Takum32}, ::Type{Takum8}) = Takum32
Base.promote_rule(::Type{Takum32}, ::Type{Takum16}) = Takum32

Base.promote_rule(::Type{Takum64}, ::Type{Takum8}) = Takum64
Base.promote_rule(::Type{Takum64}, ::Type{Takum16}) = Takum64
Base.promote_rule(::Type{Takum64}, ::Type{Takum32}) = Takum64

# IEEE 754 floating point promote rules
Base.promote_rule(::Type{Float16}, ::Type{<:AnyTakum}) = Float16
Base.promote_rule(::Type{Float32}, ::Type{<:AnyTakum}) = Float32
Base.promote_rule(::Type{Float64}, ::Type{<:AnyTakum}) = Float64

# arithmetic
Base.:(+)(x::Takum8,  y::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_addition(reinterpret(Signed, x)::Int8, reinterpret(Signed, y)::Int8)::Int8)
Base.:(+)(x::Takum16, y::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_addition(reinterpret(Signed, x)::Int16, reinterpret(Signed, y)::Int16)::Int16)
Base.:(+)(x::Takum32, y::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_addition(reinterpret(Signed, x)::Int32, reinterpret(Signed, y)::Int32)::Int32)
Base.:(+)(x::Takum64, y::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_addition(reinterpret(Signed, x)::Int64, reinterpret(Signed, y)::Int64)::Int64)

Base.:(-)(x::Takum8,  y::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_subtraction(reinterpret(Signed, x)::Int8, reinterpret(Signed, y)::Int8)::Int8)
Base.:(-)(x::Takum16, y::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_subtraction(reinterpret(Signed, x)::Int16, reinterpret(Signed, y)::Int16)::Int16)
Base.:(-)(x::Takum32, y::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_subtraction(reinterpret(Signed, x)::Int32, reinterpret(Signed, y)::Int32)::Int32)
Base.:(-)(x::Takum64, y::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_subtraction(reinterpret(Signed, x)::Int64, reinterpret(Signed, y)::Int64)::Int64)

Base.:(*)(x::Takum8,  y::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_multiplication(reinterpret(Signed, x)::Int8, reinterpret(Signed, y)::Int8)::Int8)
Base.:(*)(x::Takum16, y::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_multiplication(reinterpret(Signed, x)::Int16, reinterpret(Signed, y)::Int16)::Int16)
Base.:(*)(x::Takum32, y::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_multiplication(reinterpret(Signed, x)::Int32, reinterpret(Signed, y)::Int32)::Int32)
Base.:(*)(x::Takum64, y::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_multiplication(reinterpret(Signed, x)::Int64, reinterpret(Signed, y)::Int64)::Int64)

Base.:(/)(x::Takum8,  y::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_division(reinterpret(Signed, x)::Int8, reinterpret(Signed, y)::Int8)::Int8)
Base.:(/)(x::Takum16, y::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_division(reinterpret(Signed, x)::Int16, reinterpret(Signed, y)::Int16)::Int16)
Base.:(/)(x::Takum32, y::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_division(reinterpret(Signed, x)::Int32, reinterpret(Signed, y)::Int32)::Int32)
Base.:(/)(x::Takum64, y::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_division(reinterpret(Signed, x)::Int64, reinterpret(Signed, y)::Int64)::Int64)

Base.:(^)(x::Takum8,  y::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_power(reinterpret(Signed, x)::Int8, reinterpret(Signed, y)::Int8)::Int8)
Base.:(^)(x::Takum16, y::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_power(reinterpret(Signed, x)::Int16, reinterpret(Signed, y)::Int16)::Int16)
Base.:(^)(x::Takum32, y::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_power(reinterpret(Signed, x)::Int32, reinterpret(Signed, y)::Int32)::Int32)
Base.:(^)(x::Takum64, y::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_power(reinterpret(Signed, x)::Int64, reinterpret(Signed, y)::Int64)::Int64)

function Base.:(^)(t::AnyTakum, n::Integer)
	if n == 0
		return one(t)
	else
		nabs = abs(n)

		# get the number of bit digits
		maxpow = ceil(Int, log(2, nabs))
		t2pow = t
		nshift = nabs
		res = one(t)

		for i in 0:maxpow
			# if the i'th bit is set in n multiply res with
			# t^(2^i) to accumulate it into the product
			if nshift & 1 == 1
				res *= t2pow
			end

			t2pow *= t2pow
			nshift >>= 1
		end

		return n < 0 ? inv(res) : res
	end
end

Base.:(-)(t::AnyTakum) = Base.bitcast(typeof(t), -reinterpret(Signed, t))

# TODO ^(Takum, Takum) and ^(Takum, Integer)

Base.zero(::Type{Takum8})  = Base.bitcast(Takum8,  0x00)
Base.zero(::Type{Takum16}) = Base.bitcast(Takum16, 0x0000)
Base.zero(::Type{Takum32}) = Base.bitcast(Takum32, 0x00000000)
Base.zero(::Type{Takum64}) = Base.bitcast(Takum64, 0x0000000000000000)

Base.one(::Type{Takum8})  = Base.bitcast(Takum8,  0x40)
Base.one(::Type{Takum16}) = Base.bitcast(Takum16, 0x4000)
Base.one(::Type{Takum32}) = Base.bitcast(Takum32, 0x40000000)
Base.one(::Type{Takum64}) = Base.bitcast(Takum64, 0x4000000000000000)

Base.inv(t::Takum8)  = isnar(t) ? NaR8  : Base.bitcast(Takum8, (reinterpret(Unsigned, t) ⊻ UInt8(0x7f)) + UInt8(1))
Base.inv(t::Takum16) = isnar(t) ? NaR16 : Base.bitcast(Takum16, (reinterpret(Unsigned, t) ⊻ UInt16(0x7fff)) + UInt16(1))
Base.inv(t::Takum32) = isnar(t) ? NaR32 : Base.bitcast(Takum32, (reinterpret(Unsigned, t) ⊻ UInt32(0x7fffffff)) + UInt32(1))
Base.inv(t::Takum64) = isnar(t) ? NaR64 : Base.bitcast(Takum64, (reinterpret(Unsigned, t) ⊻ UInt64(0x7fffffffffffffff)) + UInt64(1))

# comparisons
Base.:(==)(x::T, y::T) where {T <: AnyTakum} = reinterpret(Signed, x) == reinterpret(Signed, y)
Base.:(!=)(x::T, y::T) where {T <: AnyTakum} = reinterpret(Signed, x) != reinterpret(Signed, y)
Base.:(<)(x::T, y::T) where {T <: AnyTakum} = reinterpret(Signed, x) < reinterpret(Signed, y)
Base.:(<=)(x::T, y::T) where {T <: AnyTakum} = reinterpret(Signed, x) <= reinterpret(Signed, y)

Base.isequal(x::T, y::T) where {T <: AnyTakum} = (x == y)

Base.widen(::Type{Takum8}) = Takum16
Base.widen(::Type{Takum16}) = Takum32
Base.widen(::Type{Takum32}) = Takum64

Base.widemul(x::Union{Takum8, Takum16, Takum32}, y::Union{Takum8, Takum16, Takum32}) = Basen.widen(x) * Base.widen(y)

# output
function Base.show(io::IO, t::AnyTakum)
	has_type_info = typeof(t) === Base.get(io, :typeinfo, Any)
	if isnar(t)
		Base.print(io, "NaR")
	else
		has_type_info || Base.print(io, string(typeof(t)) * "(")
		@static if VERSION ≥ v"1.7"
			@Printf.printf(IOContext(io, :typeinfo=>typeof(t)), "%.*g", max(4, 1 + Base.precision(t; base = 10)), Float64(t))
		else
			@Printf.printf(IOContext(io, :typeinfo=>typeof(t)), "%f", Float64(t))
		end
		has_type_info || Base.print(io, ")")
	end
end

Printf.tofloat(t::AnyTakum) = Float64(t)

# bitstring
Base.bitstring(t::AnyTakum) = Base.bitstring(reinterpret(Unsigned, t))

# next and previous number
Base.nextfloat(t::Takum8)  = isnar(t) ? NaR8  : Base.bitcast(Takum8,  reinterpret(Unsigned, t) + UInt8(1))
Base.nextfloat(t::Takum16) = isnar(t) ? NaR16 : Base.bitcast(Takum16, reinterpret(Unsigned, t) + UInt16(1))
Base.nextfloat(t::Takum32) = isnar(t) ? NaR32 : Base.bitcast(Takum32, reinterpret(Unsigned, t) + UInt32(1))
Base.nextfloat(t::Takum64) = isnar(t) ? NaR64 : Base.bitcast(Takum64, reinterpret(Unsigned, t) + UInt64(1))

Base.prevfloat(t::Takum8)  = isnar(t) ? NaR8  : Base.bitcast(Takum8,  reinterpret(Unsigned, t) - UInt8(1))
Base.prevfloat(t::Takum16) = isnar(t) ? NaR16 : Base.bitcast(Takum16, reinterpret(Unsigned, t) - UInt16(1))
Base.prevfloat(t::Takum32) = isnar(t) ? NaR32 : Base.bitcast(Takum32, reinterpret(Unsigned, t) - UInt32(1))
Base.prevfloat(t::Takum64) = isnar(t) ? NaR64 : Base.bitcast(Takum64, reinterpret(Unsigned, t) - UInt64(1))

# math functions
Base.abs(t::AnyTakum) = (t < 0) ? -t : t
Base.abs2(t::AnyTakum) = t * t

math_functions = [
	(:sqrt,  :square_root),
	(:cbrt,  :root, :(3::Int64)),
	(:exp,   :exp),
	(:exp2,  :(2_raised)),
	(:exp10, :(10_raised)),
	(:expm1, :exp_minus_1),
	(:log,   :ln),
	(:log2,  :lb),
	(:log10, :lg),
	(:log1p, :ln_1_plus),
	(:sin,   :sin),
	(:cos,   :cos),
	(:tan,   :tan),
	(:csc,   :csc),
	(:sec,   :sec),
	(:cot,   :cot),
	(:asin,  :arcsin),
	(:acos,  :arccos),
	(:atan,  :arctan),
	(:acsc,  :arccsc),
	(:asec,  :arcsec),
	(:acot,  :arccot),
	(:sinh,  :sinh),
	(:cosh,  :cosh),
	(:tanh,  :tanh),
	(:csch,  :csch),
	(:sech,  :sech),
	(:coth,  :coth),
	(:asinh, :arsinh),
	(:acosh, :arcosh),
	(:atanh, :artanh),
	(:acsch, :arcsch),
	(:asech, :arsech),
	(:acoth, :arcoth),
]

for (takum_type, takum_integer_type) in takum_types
	for (math_function, library_math_function, arguments...) in math_functions
		@eval begin
			Base.$math_function(t::$takum_type) = Base.bitcast($takum_type,
				@ccall libtakum.$(Symbol(lowercase(string(takum_type)), "_", library_math_function))(
				reinterpret(Signed, t)::$takum_integer_type, $(arguments...))::$takum_integer_type)
		end
	end
end

# miscellaneous
Base.bswap(t::AnyTakum) = Base.bswap_int(t)

# TODO: muladd?, rem?, mod?, random, isinf alias of isfinite?, takum?

end
