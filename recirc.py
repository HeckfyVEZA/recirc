import numpy as np
import plotly.graph_objects as go
import streamlit as st
import scipy
st.set_page_config(layout="wide")
st.header("Расчёт рециркуляции")
flow_mixture = st.number_input("Расход смеси (м³/ч)", min_value=0.0, max_value=10000000.0, value=100.0, step=1.0)
def find_m_phi(t):
    if t >= 0:
        return 611.21 * np.exp((18.678 - (t / 234.5)) * (t / (257.14 + t)))
    return 611.15 * np.exp((23.036 - (t / 333.7)) * (t / (279.82 + t)))
def find_params_by_t_and_phi(t:float, phi:float, pH=101_325):
    p_n = find_m_phi(t)
    d = 622.222 * phi * .01 * p_n / (pH - phi * .01 * p_n)
    i = (1.01 + 0.00197 * d) * t + 2.493 * d
    p_v = phi * p_n * .01
    return {'p_n': round(p_n, 2),'d': round(d, 2),'i': round(i, 2),'p_v': round(p_v, 2),'pH': round(pH, 2),'phi': round(phi, 2),'t': round(t, 2),}

atmosphere_pressure = st.number_input("Атмосферное давление (мм рт.ст.)", min_value=600.0, max_value=850.0, value=760.0, step=1.0) * 133.322
st.subheader("Параметры приточного воздуха")
cols = st.columns(2)

temperature_1 = cols[0].number_input("Температура приточного воздуха (°C)", min_value=-100.0, max_value=100.0, value=20.0, step=0.1)
relative_humidity_1= cols[1].number_input("Относительная влажность приточного воздуха (%)", min_value=0.0, max_value=100.0, value=50.0, step=0.1)
params = find_params_by_t_and_phi(temperature_1, relative_humidity_1, atmosphere_pressure)
st.subheader("Параметры вытяжного воздуха")
cols_2 = st.columns(2)
temperature_2 = cols_2[0].number_input("Температура вытяжного воздуха (°C)", min_value=-100.0, max_value=100.0, value=20.0, step=0.1)
relative_humidity_2= cols_2[1].number_input("Относительная влажность вытяжного воздуха (%)", min_value=0.0, max_value=100.0, value=50.0, step=0.1)
params_2 = find_params_by_t_and_phi(temperature_2, relative_humidity_2, atmosphere_pressure)

st.subheader("Степень рециркуляции")
first_flow_percent = st.number_input("Степень рециркуляции (%)", min_value=0.0, max_value=100.0, value=50.0, step=0.1)
mixture_temp = abs(temperature_1 - temperature_2)
mixture_temp = temperature_1 + (mixture_temp * first_flow_percent) / 100 if temperature_1 < temperature_2 else temperature_1 - (mixture_temp * first_flow_percent) / 100

mixture_d = params['d'] + abs(params_2['d'] - params['d']) * first_flow_percent / 100 if params['d'] < params_2['d'] else params['d'] - abs(params_2['d'] - params['d']) * first_flow_percent / 100

mixture_i = (1.01 + 0.00197 * mixture_d) * mixture_temp + 2.493 * mixture_d

mixture_phi = 100 * atmosphere_pressure / (((622 * find_m_phi(mixture_temp)) / mixture_d) + find_m_phi(mixture_temp))

def find_intersection(fun1, fun2, x0):
    return scipy.optimize.fsolve(lambda x : fun1(x) - fun2(x), x0)

def find_equi_i(i:float):
    d = np.array(range(1, 500, 10)) / 10
    t = (i - 2.493 * d) / (1.01 + 0.00197 * d)
    return d, t

def find_equi_phi(phi:float):
    t = np.array(range(-50, 170))
    d = []
    for ti in t:
        d += [622.222 * phi * .01 * find_params_by_t_and_phi(ti, phi, atmosphere_pressure)['p_n'] / (find_params_by_t_and_phi(ti, phi, atmosphere_pressure)['pH'] - phi * .01 * find_params_by_t_and_phi(ti, phi, atmosphere_pressure)['p_n'])]
    return d, t
condensat = 0
if mixture_phi > 100:
    mixture_phi = 100.0
    old_d = mixture_d
    f_i = scipy.interpolate.interp1d(x=find_equi_i(mixture_i)[0], y=find_equi_i(mixture_i)[1], fill_value='extrapolate')
    f_phi = scipy.interpolate.interp1d(x=find_equi_phi(mixture_phi)[0], y=find_equi_phi(mixture_phi)[1], fill_value='extrapolate')
    mixture_d = list(find_intersection(f_i, f_phi, mixture_d))
    mixture_temp = list(f_i(mixture_d))[0]
    mixture_d = mixture_d[0]
    condensat = old_d - mixture_d

fig = go.Figure()
for i in range(-35, 200, 10):
    fig.add_trace(go.Scatter(x=find_equi_i(i)[0], y=find_equi_i(i)[1], mode='lines', name=str(i), line=dict(color='black', width=.3)), )
for phi in range(0, 101, 10):
    fig.add_trace(go.Scatter(x=find_equi_phi(phi)[0], y=find_equi_phi(phi)[1], mode='lines', name=str(phi), line=dict(color='black', width=.3)))

fig.add_trace(go.Scatter(x=[params['d']], y=[params['t']], mode='markers', name=f'{params["d"]}, {params["t"]}', line=dict(color='green', width=5)))
fig.add_trace(go.Scatter(x=find_equi_i(params['i'])[0], y=find_equi_i(params['i'])[1], mode='lines', name="i = " + str(params['i']), line=dict(color='green', width=.35, dash='dash')), )
fig.add_trace(go.Scatter(x=find_equi_phi(params['phi'])[0], y=find_equi_phi(params['phi'])[1], mode='lines', name="φ = " + str(params['phi']), line=dict(color='green', width=.35, dash='dash')))
fig.add_trace(go.Scatter(x=[params['d'] for l in range(0, 200, 10)], y=[l-50 for l in range(0, 200, 10)], mode='lines', name="d = " + str(params['d']), line=dict(color='green', width=.35, dash='dash')))
fig.add_trace(go.Scatter(x=[l for l in range(0, 200, 10)], y=[params['t'] for l in range(0, 200, 10)], mode='lines', name= "t = " + str(params['t']), line=dict(color='green', width=.35, dash='dash')))


fig.add_trace(go.Scatter(x=[params_2['d']], y=[params_2['t']], mode='markers', name=f'{params_2["d"]}, {params_2["t"]}', line=dict(color='blue', width=5, dash='dash')))
fig.add_trace(go.Scatter(x=find_equi_i(params_2['i'])[0], y=find_equi_i(params_2['i'])[1], mode='lines', name="i = " + str(params_2['i']), line=dict(color='blue', width=.35, dash='dash')), )
fig.add_trace(go.Scatter(x=find_equi_phi(params_2['phi'])[0], y=find_equi_phi(params_2['phi'])[1], mode='lines', name="φ = " + str(params_2['phi']), line=dict(color='blue', width=.35)))
fig.add_trace(go.Scatter(x=[params_2['d'] for l in range(0, 200, 10)], y=[l-50 for l in range(0, 200, 10)], mode='lines', name="d = " + str(params_2['d']), line=dict(color='blue', width=.35, dash='dash')))
fig.add_trace(go.Scatter(x=[l for l in range(0, 200, 10)], y=[params_2['t'] for l in range(0, 200, 10)], mode='lines', name= "t = " + str(params_2['t']), line=dict(color='blue', width=.35, dash='dash')))

fig.add_trace(go.Scatter(x=[params['d'], params_2['d']], y=[params['t'], params_2['t']], mode='lines', name=f'{params["d"]}, {params["t"]}', line=dict(color='black', width=1, dash='dash')))

fig.add_trace(go.Scatter(x=[mixture_d], y=[mixture_temp], mode='markers', name=f'{mixture_d}, {mixture_temp}', line=dict(color='red', width=5)))
fig.add_trace(go.Scatter(x=find_equi_i(mixture_i)[0], y=find_equi_i(mixture_i)[1], mode='lines', name="i = " + str(mixture_i), line=dict(color='red', width=1, dash='dash')), )
fig.add_trace(go.Scatter(x=find_equi_phi(mixture_phi)[0], y=find_equi_phi(mixture_phi)[1], mode='lines', name="φ = " + str(mixture_phi), line=dict(color='red', width=1, dash='dash')))
fig.add_trace(go.Scatter(x=[mixture_d for l in range(0, 200, 10)], y=[l-50 for l in range(0, 200, 10)], mode='lines', name="d = " + str(mixture_d), line=dict(color='red', width=1, dash='dash')))
fig.add_trace(go.Scatter(x=[l for l in range(0, 200, 10)], y=[mixture_temp for l in range(0, 200, 10)], mode='lines', name= "t = " + str(mixture_temp), line=dict(color='red', width=1, dash='dash')))

fig.update_layout(title='i-d diagramm', xaxis_title='d, г/м³', yaxis_title='t, °С', xaxis_range=[0, 45], yaxis_range=[-35, 80], showlegend=False, width=1000, height=1000, xaxis=dict(showgrid=True, zeroline=True), yaxis=dict(showgrid=True, zeroline=True))

figcols = st.columns(3)
figcols[0].plotly_chart(fig)
with figcols[2]:
    st.write(f"Относительная влажность смеси: {round(mixture_phi, 2)} %")
    st.write(f"Влагосодержание смеси: {round(mixture_d, 2)} г/м³")
    st.write(f"Энтальпия смеси: {round(mixture_i, 2)} кДж/кг")
    st.write(f"Температура смеси: {round(mixture_temp, 2)} °С")
    st.write(f"Выпадет конденсата: {round(condensat * flow_mixture, 2)} г/ч")
