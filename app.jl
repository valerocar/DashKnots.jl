using Pkg
Pkg.activate(".")
using PlotlyJS, Zygote, LinearAlgebra
using Dash, DashCoreComponents, DashHtmlComponents
using DashBootstrapComponents

function frenet_frame(x, y, z)
    ẋ(t) = gradient(x, t)[1]
    ẏ(t) = gradient(y, t)[1]
    ż(t) = gradient(z, t)[1]

    ẍ(t) = hessian(x, t)
    ÿ(t) = hessian(y, t)
    z̈(t) = hessian(z, t)

    α̇(t) = [ẋ(t), ẏ(t), ż(t)]
    α̈(t) = [ẍ(t), ÿ(t), z̈(t)]
    function out(t)
        T = α̇(t)
        T = T / norm(T)
        ad = α̈(t)
        N = ad - dot(ad, T)T
        N = N / norm(N)
        B = cross(T, N)
        T, N, B
    end
end

function uv_surface(uv_domain, ϕ; uv_res = [40, 40])
    us = LinRange(uv_domain[1], uv_domain[2], uv_res[1])
    vs = LinRange(uv_domain[3], uv_domain[4], uv_res[2])
    u_res, v_res = uv_res
    xs = zeros(u_res, v_res)
    ys = zeros(u_res, v_res)
    zs = zeros(u_res, v_res)
    for (i, u) in enumerate(us)
        for (j, v) in enumerate(vs)
            xs[i, j], ys[i, j], zs[i, j] = ϕ(u, v)
        end
    end
    surface(x = xs, y = ys, z = zs, showscale = false)
end

function xyz(ps)
    x = [p[1] for p in ps]
    y = [p[2] for p in ps]
    z = [p[3] for p in ps]
    return x, y, z
end

function torus_knot(m::Int64, n::Int64; R = 2, r = 1, uv_res = [220, 220])
    ϕ(u, v) = [(R + r * cos(u)) * cos(v), (R + r * cos(u)) * sin(v), r * sin(u)]
    α(t) = ϕ(m * t, n * t)
    x(t) = α(t)[1]
    y(t) = α(t)[2]
    z(t) = α(t)[3]

    F = frenet_frame(x, y, z)
    function ψ(u, v)
        T, N, B = F(u)
        ρ = r / 3
        α(u) + ρ * cos(v) * N + ρ * sin(v) * B
    end
    uv_surface([-π, 1.02π, -π, π], ψ, uv_res = uv_res)
end
description = dcc_markdown(raw"Torus knots are obtained by winding a curve around a torus 
 ($m$ times in parallel direction and $n$ times in the meridian direction). More concretely, if

 $$\varphi(u,v) = ((R+r\cos(u))\cos(v), (R+r\cos(u))\sin(v), r\sin(u))$$

 is the torus parametrization, the curve

 $$\alpha(t) = \varphi(mt,nt)\text{ for }0\leq t \leq 2\pi $$ 

is said to be knot of the type $(m,n)$. These knots can be made into surfaces by thickening them using the Frenet frame over $\alpha$
 ")
sd = 4.0
scene3d = Dict(
    :xaxis => Dict(:range => [-sd, sd], :autorange => false),
    :yaxis => Dict(:range => [-sd, sd], :autorange => false),
    :zaxis => Dict(:range => [-sd, sd], :autorange => false),
    :aspectratio => Dict(:x => 1, :y => 1, :z => 1),
)

layout = Layout(scene = scene3d, width = 600, height = 600)

mathjax = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML"
app = dash(external_scripts = [mathjax], external_stylesheets = [dbc_themes.BOOTSTRAP])
app.title = "Dash Torus Knots"

knot = torus_knot(2, 3)
graph = dbc_spinner(dcc_graph(id = "graph", figure = plot(knot, layout)))

slider_m = dcc_slider(
    id = "slider_m",
    min = 1,
    max = 7,
    step = 1,
    value = 3,
    marks = Dict([i => ("m=$(i)") for i = 1:10]),
)
slider_n = dcc_slider(
    id = "slider_n",
    min = 1,
    max = 7,
    step = 1,
    value = 2,
    marks = Dict([i => ("n=$(i)") for i = 1:10]),
)



input = dbc_row([dbc_col(slider_m, width = Dict("size" => 3, "order" => "first", "offset" => 3)),
 dbc_col(slider_n, width = Dict("size" => 3, "order" => "last", "offset" => 0))])
app.layout = dbc_container(
    [
        html_h1(html_center("Torus Knots"))
        html_br()
        description
        html_center(graph)
        html_center(input)
        html_br()
        html_br()
    ],
)

callback!(
    app,
    Output("graph", "figure"),
    Input("slider_m", "value"),
    Input("slider_n", "value"),
) do m, n
    knot = torus_knot(m, n)
    return plot(knot, layout)
end
run_server(app, "0.0.0.0", debug = true)
