-- Hyperbolic functions not available in Lua 5.3 anymore, implement them
-- on our own here (naive implementation, use Fortran instead for more accurate
-- implementations).

function sinh(x)
  return((math.exp(2*x)-1)/(2*math.exp(x)))
end

function cosh(x)
  return((math.exp(2*x)+1)/(2*math.exp(x)))
end

function tanh(x)
  return((math.exp(2*x)-1)/(math.exp(2*x)+1))
end
