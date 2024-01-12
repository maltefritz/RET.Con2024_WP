# -*- coding: utf-8 -*-
"""Functions for calculation of economical and ecological parameters.

Created on Mon Feb  1 11:43:39 2021

@author: Jonas Freißmann
"""

import numpy as np
import pandas as pd

try:
    from scipy.optimize import curve_fit
except ImportError:
    print(
        'WARNING: function "curve_fit" not available, because scipy could not '
        + 'be imported'
        )

# import matplotlib.pyplot as plt

def calc_bwsf(i, n):
    """Berechne Barwert Summenfaktor.
    
    Parameters:
    -----------
    i : float
        Kapitalzins als rationale Zahl (nicht Prozent)

    n : int
        Lebensdauer des Investments in Jahren
    """
    q = 1+i
    return (q**n - 1)/(q**n * (q - 1))


def npv(invest, cashflow, i=0.05, n=20):
    """Konstantin 2013, Markus [29].

    npv:        Kapitalwert (netpresentvalue)
    invest:     Investitionsausgaben zum Zeitpunkt t=0
    cashflow:   Differenz aller Einnahmen und Ausgaben (Zahlungsströme)
    i:          Kalkulationszinssatz
    n:          Betrachtungsdauer
    bwsf:       Barwert Summenfaktor
    """
    q = 1+i
    bwsf = (q**n - 1)/(q**n * (q - 1))

    npv = -invest + bwsf * cashflow
    return npv


def LCOH_old(invest, cashflow, Q, i=0.05, n=20):
    """Konstantin 2013, Markus [29].

    LCOH        Wärmegestehungskosten
    invest:     Investitionsausgaben zum Zeitpunkt t=0
    bwsf:       Barwert Summenfaktor
    cashflow:   Differenz aller Einnahmen und Ausgaben (Zahlungsströme)
                innerhalb des betrachteten Jahres
    Q:          Gesamte bereitgestellte Wärmemenge pro Jahr
    i:          Kalkulationszinssatz
    n:          Betrachtungsdauer
    """
    q = 1 + i
    bwsf = (q**n - 1)/(q**n * (q - 1))

    LCOH = abs(-invest * bwsf**(-1) + cashflow) / Q
    return LCOH


def LCOH(invest, cost, Q, revenue=0, i=0.05, n=20):
    """Konstantin 2013, Markus [29].

    LCOH        Wärmegestehungskosten
    invest:     Investitionsausgaben zum Zeitpunkt t=0
    bwsf:       Barwert Summenfaktor
    cashflow:   Differenz aller Einnahmen und Ausgaben (Zahlungsströme)
                innerhalb des betrachteten Jahres
    Q:          Gesamte bereitgestellte Wärmemenge pro Jahr
    i:          Kalkulationszinssatz
    n:          Betrachtungsdauer
    """
    q = 1 + i
    bwsf = (q**n - 1)/(q**n * (q - 1))

    LCOH = (invest + bwsf * (cost - revenue))/(bwsf * Q)
    return LCOH


def emission_calc(data_emission):
    """Calculate the emissions compared with overall and displacement mix."""
    co2 = pd.read_csv('input\\emissions2016.csv', sep=";")
    co2['Date'] = pd.to_datetime(co2['Date'], format='%d.%m.%Y %H:%M')
    co2.set_index('Date', inplace=True)

    Em_om = list()
    Em_dm = list()

    e_fuel = 0.2012    # in t/MWh aus "Emissionsbewertung Daten"

    for idx_em, idx_co2 in zip(data_emission.index, co2.index):
        Em_om.append(
            data_emission.loc[idx_em, 'H_source'] * e_fuel
            + data_emission.loc[idx_em, 'P_source'] * co2.loc[idx_co2, 'EF om']
            - (data_emission.loc[idx_em, 'P_spot_market']
               * co2.loc[idx_co2, 'EF om'])
            )

        Em_dm.append(
            data_emission.loc[idx_em, 'H_source'] * e_fuel
            + data_emission.loc[idx_em, 'P_source'] * co2.loc[idx_co2, 'EF dm']
            - (data_emission.loc[idx_em, 'P_spotmarket']
               * co2.loc[idx_co2, 'EF dm'])
            )

    dfEm = pd.DataFrame({'date': data_emission.index,
                         'overall mix': Em_om,
                         'displacement mix': Em_dm})
    dfEm.set_index('date', inplace=True)

    return dfEm


def emission_calc(data_all, data, param):
    """
    Calculate the emissions compared with overall and displacement mix.

    WARNING! Only works with the primary network for now.

    Parameters
    ----------
    data_all : pandas.DataFrame
        result DataFrame containing all flows for each time step.

    data : pandas.DataFrame
        csv file of user defined time dependent parameters containing the
        specific emission factors.

    param : dict
        JSON parameter file of user defined constants.
    """
    print('WARNING! Only works with the primary network for now.')
    # data_all['Emissions OM'] = [0 for _ in range(len(data_all.index))]
    # data_all['Emissions DM'] = [0 for _ in range(len(data_all.index))]
    data_all['Emissions OM'] = 0
    data_all['Emissions DM'] = 0

    if 'P_sub_source' in data_all.columns:
        data_all['Sub Emissions OM'] = 0
        data_all['Sub Emissions DM'] = 0

    if 'H_source' in data_all.columns:
        data_all['Emissions OM'] += (
            data_all['H_source'] * param['param']['ef_gas']
            )
        data_all['Emissions DM'] += (
            data_all['H_source'] * param['param']['ef_gas']
            )

    if 'H_bio_source' in data_all.columns:
        data_all['Emissions OM'] += (
            data_all['H_bio_source'] * param['param']['ef_biogas']
            )
        data_all['Emissions DM'] += (
            data_all['H_bio_source'] * param['param']['ef_biogas']
            )

    if 'P_source' in data_all.columns:
        data_all['Emissions OM'] += data_all['P_source'] * data['ef_om']
        data_all['Emissions DM'] += data_all['P_source'] * data['ef_dm']
    if 'P_sub_source' in data_all.columns:
        data_all['Sub Emissions OM'] += data_all['P_sub_source'] * data['ef_om']
        data_all['Sub Emissions DM'] += data_all['P_sub_source'] * data['ef_dm']


    if 'P_spotmarket' in data_all.columns:
        data_all['Emissions OM'] -= data_all['P_spotmarket'] * data['ef_om']
        data_all['Emissions DM'] -= data_all['P_spotmarket'] * data['ef_dm']

    return data_all


def invest_sol(A, col_type=''):
    """Pehnt et al. 2017, Markus [38].

    A:                Kollektorfläche der Solarthermie
    col_type:         Kollektortyp der Solarthermie
    specific_coasts:  Spezifische Kosten
    invest:           Investitionskosten
    """
    if col_type == 'flat':
        specific_costs = -34.06 * np.log(A) + 592.48
        invest = A * specific_costs
        return invest
    elif col_type == 'vacuum':
        specific_costs = -40.63 * np.log(A) + 726.64
        invest = A * specific_costs
        return invest
    else:
        raise ValueError(
            "Choose a valid collector type: 'flat' or 'vacuum'"
            )


def invest_stes(Q):
    """Investment calculation for seasonal thermal energy storages.

    Q:              Kapazität des Speichers in MWh
    sponsorship:    Förderung des Speichers (durch Bundesamt für
                    Wirtschaft und Ausfuhrkontrolle [10])
    q_V:            spez. volumetrische Energie
    """
    # Kostendegression STES
    def potential_func(x, a, b):
        return a * x ** b

    if Q == 0:
        return 0

    x = [500, 5000, 62000]
    y = [320, 110, 2359594/62000]

    params, params_covariance = curve_fit(potential_func, x, y)

    # V = Q / 0.3
    q_V = 0.07    # MWh/m³, unsere Annahme (siehe Whiteboardbild)
    V_stes = Q / q_V
    Q_specific_costs = params[0] * V_stes ** params[1]
    stes_invest = V_stes * Q_specific_costs

    if V_stes > 50:
        sponsorship = 250 * V_stes
        if sponsorship > 10e6:
            sponsorship = 10e6
        if sponsorship > 0.3 * stes_invest:
            sponsorship = 0.3 * stes_invest
        stes_invest -= sponsorship

    return stes_invest


def chp_bonus(P, use_case):
    """Calculate chp bonus based on nominal power output in.

    P:           nomimal power output of chp unit in kW (int or float)
    use_case:    either 'grid' or 'self-sufficient' (str)
    bonus:       calculated chp bonus in ct/kWh (float)
    """
    if P == 0:
        return 0
    if use_case == 'grid':
        bonus_intervals = [8.0, 6.0, 5.0, 4.4, 3.4]
    elif use_case == 'self-sufficient':
        bonus_intervals = [4.0, 3.0, 2.0, 1.5, 1.0]
    else:
        print('No valid use case given.')

    # Defined intervals for power output (KWKG)
    P_intervals = [0.0, 50.0, 100.0, 250.0, 2000.0]

    # Calculate steps between power output intervals
    P_steps = list()
    for i in range(1, len(P_intervals)):
        P_steps += [P_intervals[i] - P_intervals[i-1]]

    # Check at which index the last power value is bigger than P
    idx = len(P_intervals) - 1
    for i in range(len(P_intervals)):
        if P < P_intervals[i]:
            idx = i - 1
            break

    # Add the weighted bonus values for all complete intervals
    bonus_weighted = 0
    for i in range(idx):
        bonus_weighted += P_steps[i] * bonus_intervals[i]

    # Add the weighted bonus for the last incomplete interval
    bonus_weighted += (P - P_intervals[idx]) * bonus_intervals[idx]

    # Calculate the nominal bonus by deviding by the sum of weights (P)
    bonus = bonus_weighted / P

    return bonus


def bew_op_bonus(SCOP, el_source_type='conventional'):
    """
    Calculate operating cost BEW bonus for heat pumps.

    Result in €/MWh and capped at 92€/MWh. Bonus shall not be greater than 90%
    of real electricity cost. SCOP has to be greater than or equal to 2.5 to
    qualify 

    Source of formula:
    https://www.bafa.de/SharedDocs/Downloads/DE/Energie/bew_merkblatt_antragstellung_m4.pdf?__blob=publicationFile&v=2

    Parameters
    ----------

    SCOP : float
        Seasonal Coefficient of Performance of the heat pump.

    el_source_type : str
        Type of electricity source. Either 'conventional' for grid or
        'renewable' for direct use of renewable electricity sources.
        Defaults to 'conventional'
    """
    if el_source_type == 'conventional':
        bonus = (5.5 - (6.8 - 17/SCOP) * 0.75) * SCOP/(SCOP - 1)
    elif el_source_type == 'renewable':
        bonus = 3 - (8/2.5 - 8/SCOP) * 0.75

    return min(bonus*10, 92.00)


if __name__ == '__main__':
    # %% Invest cost linearization of solar thermal collectors
    # plot_cost = True
    # A_range = np.linspace(1, 100000, 500)
    # I_range = invest_sol(A_range, col_type='flat') / 1e6

    # polreg = np.poly1d(np.polyfit(A_range, I_range, 1))
    # I_lin_range = polreg(A_range)

    # if plot_cost:
    #     import matplotlib.pyplot as plt
    #     fig, ax = plt.subplots(figsize=(8, 6))

    #     ax.scatter(A_range, I_range, lw=0.05)
    #     ax.plot(A_range, I_lin_range, color='r')

    #     ax.set_xlabel('Kollektorfläche in m²')
    #     ax.set_ylabel('Investitionskosten in Mio. €')

    #     ax.grid()

    #     plt.show()

    # c0 = round(polreg.coefficients[1]*1e6, 2)
    # c1 = round(polreg.coefficients[0]*1e6, 2)
    # print(f'c0 = {c0}, c1 = {c1}')

    # scop = 3.0
    # print(bew_op_bonus(scop, el_source_type='conventional'))
    # print(bew_op_bonus(scop, el_source_type='renewable'))

    # Q_sttes_range = np.linspace(6*24, 353*24, 10)
    # for Q in Q_sttes_range:
    #     print(round(invest_stes(Q), 2))
    # m = (
    #     (invest_stes(Q_sttes_range[-1]) - invest_stes(Q_sttes_range[0]))
    #     / (Q_sttes_range[-1] - Q_sttes_range[0])
    #     )
    # b = invest_stes(Q_sttes_range[0]) - m * Q_sttes_range[0]
    # print(f'm: {m:.2f} €/MWh, b: {b:.2f} €')
    # result: m: 236.49 €/MWh, b: 207314.97 €

    Q_stes_range = np.linspace(100, 100000, 10)
    for Q in Q_stes_range:
        print(round(invest_stes(Q), 2))
    m = (
        (invest_stes(Q_stes_range[-1]) - invest_stes(Q_stes_range[0]))
        / (Q_stes_range[-1] - Q_stes_range[0])
        )
    b = invest_stes(Q_stes_range[0]) - m * Q_stes_range[0]
    print(f'm: {m:.2f} €/MWh, b: {b:.2f} €')
    # result: m: 236.49 €/MWh, b: 207314.97 €
#     invest = [
#         7661674.2799370885, 7656529.477196013, 7684854.4973669015
#     ]

#     cost = [
#         1460292.0152750215, 1412236.673756644, 1386007.325025606
#     ]

#     Q_total = [21781.65304, 21781.65304, 21781.65304]

#     revenues = [
#         1169584.0867732146, 1131152.5787058724, 1110914.5070661784
#     ]

#     lcoh_res = []
#     for inv, co, Q, rev in zip(invest, cost, Q_total, revenues):
#         lcoh_res += [LCOH(inv, co, Q, rev)]

#     print(lcoh_res)

#     import matplotlib.pyplot as plt
#     import SWSHplotting as shplt

#     shplt.init_params()

#     fig, ax = plt.subplots(figsize=(12, 8))

#     ax.bar([0, 1, 2], lcoh_res, color=shplt.znes_colors()['darkblue'])

#     ax.grid(':', axis='y')

#     ax.set_ylabel('Wärmegestehungskosten $LCOH$ in $€/MWh$')
#     ax.set_xticks([0, 1, 2])
#     ax.set_xticklabels(['Einfacher Kreis - R1234yf', 'Economizer - R1234yf', 'Economizer - Ammoniak'])

#     for i, lcoh in enumerate(lcoh_res):
#         ax.annotate(
#             f'{lcoh:.2f}', (i, lcoh * 1.05), ha='center', va='center', color='k'
#         )

#     ax.set_ylim(top=max(lcoh_res)*1.1)

#     # plt.savefig('lcoh_vergleich.png')
#     plt.show()

    # import matplotlib.pyplot as plt

    # znesstyle = True

    # if znesstyle:
    #     import SWSHplotting as shplt
    #     shplt.init_params()

    # chpboni = [
    #     chp_bonus(P_N, 'grid')*1000/100  for P_N in np.linspace(1e2, 3.5e5, 500)
    #     ]

    # print(f'KWK-Bonus @ 50MW: {chp_bonus(0.5e6, "grid")*1000/100:.2f} €/MWh')
    # print(f'KWK-Bonus @ 100MW: {chp_bonus(1e6, "grid")*1000/100:.2f} €/MWh')
    # print(f'KWK-Bonus @ 200MW: {chp_bonus(2e6, "grid")*1000/100:.2f} €/MWh')
    # print(f'KWK-Bonus @ 350MW: {chp_bonus(3.5e6, "grid")*1000/100:.2f} €/MWh')

    # fig, ax = plt.subplots(figsize=(8, 5))

    # if znesstyle:
    #     ax.axhline(
    #         3.4*1000/100, color=shplt.znes_colors()['red'], lw=2, ls='--',
    #         label='Letztes Bonusintervall (34 $€/MWh$)'
    #         )
    #     ax.plot(
    #         np.linspace(0.1, 300, 500), chpboni, lw=2,
    #         color=shplt.znes_colors()['darkblue']
    #         )
    #     ax.legend(fontsize=14)
    # else:
    #     ax.axhline(3.4*1000/100, color='tab:red', lw=2, ls='--')
    #     ax.plot(np.linspace(0.1, 300, 500), chpboni, lw=2)

    # ax.grid()

    # ax.set_ylim(bottom=0)

    # ax.set_xlabel('Elektrische Nennleistung in $MW$')
    # ax.set_ylabel('KWK-Bonus in $€/MWh$')

    # plt.tight_layout()
    # # if znesstyle:
    # #     plt.savefig('KWK_Bonus_über_Nennleistung.pdf')

    # plt.show()