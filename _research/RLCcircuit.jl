using DifferentialEquations
using Plots

# Function that takes the current state of the system and returns the derivative of the state
function rlc_circuit(du, u, p, t)
    Q, dQ = u         # Q: charge on the capacitor
    L, R, C, V = p
    
    du[1] = dQ
    du[2] = (V - Q/C - R*dQ)/L  # The differential equation for RLC circuit is:
                                # L(d^2Q/dt^2) + R(dQ/dt) + Q/C = V(t)
        # du: the derivative of the state (output)
end

Q0 = 0.0
dQ0 = 0.0
u0 = [Q0, dQ0]

L = 1.0     # L: inductance
R = 0.5     # R: resistance
C = 0.1     # C: capacitance
V = 1.0
p = (L, R, C, V)

tspan = (0.0, 10.0)


prob = ODEProblem(rlc_circuit, u0, tspan, p)
sol = solve(prob)

"""
plot(sol.t, [sol.u[i][1] for i in 1:length(sol.u)], label="Q(t)")
plot!(sol.t, [sol.u[i][2] for i in 1:length(sol.u)], label="dQ/dt(t)")
xlabel!("Time")
ylabel!("Charge")
"""

plot(sol.t, [sol.u[i][1] for i in 1:length(sol.u)], label="Charge on Capacitor")
xlabel!("Time")
ylabel!("Charge")