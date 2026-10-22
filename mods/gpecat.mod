TITLE T-type calcium current (Terman et al. 2002)
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
}
NEURON {
    SUFFIX gpecat
    USEION ca READ eca WRITE ica
    RANGE ggpecatbar, ica
}
PARAMETER { ggpecatbar = 0.0005 (S/cm2) }
STATE { m h }
ASSIGNED {
    v (mV)
    eca (mV)
    ica (mA/cm2)
    m_inf
    h_inf
    tau_m (ms)
    tau_h (ms)
}
BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = ggpecatbar*m*m*h*(v - eca)
}
DERIVATIVE states {
    rates(v)
    m' = (m_inf - m)/tau_m
    h' = (h_inf - h)/tau_h
}
INITIAL {
    rates(v)
    m = m_inf
    h = h_inf
}
PROCEDURE rates(v(mV)) {
    m_inf = 1/(1+exp(-(v+63)/7.8))
    h_inf = 1/(1+exp((v+80)/6.6))
    tau_m = 1.5 + 25/(1+exp((v+60)/5))
    tau_h = 10 + 50/(1+exp((v+70)/5))
}
