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

# worst case precision n-12 (TODO: instance precision)
function _precision(x::AnyTakum; base::Integer)
	base > 1 || throw(DomainError(base, "`base` cannot be less than 2."))
	p = 0 # TODO: Takum8:0, Takum16:4, Takum32:20, Takum64:52
	return base == 2 ? Int(p) : floor(Int, p / log2(base))
end
Base.precision(::Type{AnyTakum}; base::Integer=2) = _precision(T, base)
Base.precision(t::AnyTakum; base::Integer=2) = precision(typeof(t); base)

# rounding
Base.round(t::AnyTakum) = typeof(t)(Base.round(Float64(t)))

Base.round(t::AnyTakum, r::RoundingMode{:ToZero})  = typeof(t)(Base.trunc(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Down})    = typeof(t)(Base.floor(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Up})      = typeof(t)(Base.ceil(Float64(t)))
Base.round(t::AnyTakum, r::RoundingMode{:Nearest}) = typeof(t)(Base.round(Float64(t)))

Base.trunc(t::AnyTakum) = Base.signbit(t) ? Base.ceil(t) : Base.floor(t)

# type limits
Base.typemin(::Type{Takum8})  = NaR8
Base.typemin(::Type{Takum16}) = NaR16
Base.typemin(::Type{Takum32}) = NaR32
Base.typemin(::Type{Takum64}) = NaR64

Base.floatmin(::Type{Takum8})  = Base.bitcast(Takum8,  0x81)
Base.floatmin(::Type{Takum16}) = Base.bitcast(Takum16, 0x8001)
Base.floatmin(::Type{Takum32}) = Base.bitcast(Takum32, 0x80000001)
Base.floatmin(::Type{Takum64}) = Base.bitcast(Takum64, 0x8000000000000001)

Base.typemax(::Type{Takum8})  = Base.bitcast(Takum8,  0x7f)
Base.typemax(::Type{Takum16}) = Base.bitcast(Takum16, 0x7fff)
Base.typemax(::Type{Takum32}) = Base.bitcast(Takum32, 0x7fffffff)
Base.typemax(::Type{Takum64}) = Base.bitcast(Takum64, 0x7fffffffffffffff)

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
for takum_type in (Takum8, Takum16, Takum32, Takum64)
	for int_type in (Bool, UInt8, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64)
		@eval begin
			(::Type{$takum_type})(i::($int_type)) = $takum_type(Float64(i))
			Base.promote_rule(::Type{$takum_type}, ::Type{$int_type}) = $takum_type
		end
	end
end

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
Takum8(i::Integer)  = Base.convert(Takum8,  Base.convert(Float64, i))
Takum16(i::Integer) = Base.convert(Takum16, Base.convert(Float64, i))
Takum32(i::Integer) = Base.convert(Takum32, Base.convert(Float64, i))
Takum64(i::Integer) = Base.convert(Takum64, Base.convert(Float64, i))

Base.unsafe_trunc(T::Type{<:Integer}, t::Takum8)  = Base.unsafe_trunc(T, Float64(t))
Base.unsafe_trunc(T::Type{<:Integer}, t::Takum16) = Base.unsafe_trunc(T, Float64(t))
Base.unsafe_trunc(T::Type{<:Integer}, t::Takum32) = Base.unsafe_trunc(T, Float64(t))
Base.unsafe_trunc(T::Type{<:Integer}, t::Takum64) = Base.unsafe_trunc(T, Float64(t))

for int_type in (Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128)
	for takum_type in (Takum8, Takum16, Takum32, Takum64)
		@eval begin
			function Base.round(::Type{$int_type}, t::$takum_type)
				if x == -1
					return -Base.one($takum_type)
				elseif x == 0
					return Base.zero($takum_type)
				elseif x == 1
					return Base.one($takum_type)
				else
					throw(InexactError(:round, $int_type, t))
				end
			end
		end
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
function Base.show(io::IO, t::Takum8)
	has_type_info = Takum8 === Base.get(io, :typeinfo, Any)
	if isnar(t)
		Base.print(io, "NaR")
	else
		has_type_info || Base.print(io, "Takum8(")
		Base.show(IOContext(io, :typeinfo=>Float16), Float16(t))
		has_type_info || Base.print(io, ")")
	end
end
function Base.show(io::IO, t::Takum16)
	has_type_info = Takum16 === Base.get(io, :typeinfo, Any)
	if isnar(t)
		Base.print(io, "NaR")
	else
		has_type_info || Base.print(io, "Takum16(")
		Base.show(IOContext(io, :typeinfo=>Float32), Float32(t))
		has_type_info || Base.print(io, ")")
	end
end
function Base.show(io::IO, t::Takum32)
	has_type_info = Takum32 === Base.get(io, :typeinfo, Any)
	if isnar(t)
		Base.print(io, "NaR")
	else
		has_type_info || Base.print(io, "Takum32(")
		Base.show(IOContext(io, :typeinfo=>Float64), Float64(t))
		has_type_info || Base.print(io, ")")
	end
end
function Base.show(io::IO, t::Takum64)
	has_type_info = Takum64 === Base.get(io, :typeinfo, Any)
	if isnar(t)
		Base.print(io, "NaR")
	else
		has_type_info || Base.print(io, "Takum64(")
		Base.show(IOContext(io, :typeinfo=>Float64), Float64(t))
		has_type_info || Base.print(io, ")")
	end
end

Printf.tofloat(t::Takum8)  = Float16(t)
Printf.tofloat(t::Takum16) = Float32(t)
Printf.tofloat(t::Takum32) = Float64(t)
Printf.tofloat(t::Takum64) = Float64(t)

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

# TODO math functions
Base.abs(t::AnyTakum) = (t < 0) ? -t : t
Base.abs2(t::AnyTakum) = t * t

Base.sqrt(t::Takum8)  = Base.bitcast(Takum8,  @ccall libtakum.takum8_square_root(reinterpret(Signed, t)::Int8)::Int8)
Base.sqrt(t::Takum16) = Base.bitcast(Takum16, @ccall libtakum.takum16_square_root(reinterpret(Signed, t)::Int16)::Int16)
Base.sqrt(t::Takum32) = Base.bitcast(Takum32, @ccall libtakum.takum32_square_root(reinterpret(Signed, t)::Int32)::Int32)
Base.sqrt(t::Takum64) = Base.bitcast(Takum64, @ccall libtakum.takum64_square_root(reinterpret(Signed, t)::Int64)::Int64)

for func in (:cbrt, :exp, :exp2, :exp10, :expm1,
          :log, :log2, :log10, :log1p,
          :sin, :cos, :tan, :csc, :sec, :cot,
          :asin, :acos, :atan, :acsc, :asec, :acot,
          :sinh, :cosh, :tanh, :csch, :sech, :coth,
          :asinh, :acosh, :atanh, :acsch, :asech, :acoth)
  @eval begin
     # TODO call into library
     Base.$func(t::Takum8)  = Takum8(Base.$func(Float16(t)))
     Base.$func(t::Takum16) = Takum16(Base.$func(Float32(t)))
     Base.$func(t::Takum32) = Takum32(Base.$func(Float64(t)))
     Base.$func(t::Takum64) = Takum64(Base.$func(Float64(t)))
  end
end

# miscellaneous
Base.bswap(t::AnyTakum) = Base.bswap_int(t)

# TODO: muladd?, rem?, mod?, random, isinf alias of isfinite?, takum?

end
