import numpy as np
import vegas as vg
import pandas as pd

from DarkNews.phase_space import three_body_decay
from DarkNews.MC import run_vegas
from DarkNews.MC import get_samples
from DarkNews import const

from nuflux import integrands


class MC_events(object):
    def __init__(
        self,
        model="LOmudecay_unpol",
        Mparent=const.m_mu,
        Mdaughter=const.m_e,
        pmin=0.010,
        pmax=10.0,
        tmin=0.0,  # rad
        tmax=1.0,  # rad
        include_beamdiv=True,
        helicity=-1,
        mnu1=0,
        mnu2=0,
        beam_p0=3.8,  # GeV
        beam_dpop=0.10,  # %
        beam_theta0=0.0,  # rad
        beam_dtheta=0.005,  # rad
        NINT=10,
        NINT_warmup=10,
        NEVAL=1e5,
        NEVAL_warmup=1e4,
    ):
        # set target properties
        self.pmin = pmin
        self.pmax = pmax
        self.tmin = tmin
        self.tmax = tmax
        self.Mparent = Mparent
        self.Mdaughter = Mdaughter
        self.mnu1 = mnu1
        self.mnu2 = mnu2
        self.helicity = helicity
        self.model = model
        self.include_beamdiv = include_beamdiv

        self.beam_p0 = beam_p0
        self.beam_dpop = beam_dpop
        self.beam_theta0 = beam_theta0
        self.beam_dtheta = beam_dtheta

        self.NINT = NINT
        self.NINT_warmup = NINT_warmup
        self.NEVAL = NEVAL
        self.NEVAL_warmup = NEVAL_warmup

    def get_MC_events(self):
        if self.model == "LOmudecay_unpol":
            # BATCH SAMPLE INTEGRAN OF INTEREST
            DIM = 6
            batch_f = integrands.LOmudecay_unpol(dim=DIM, MC_case=self)
            integ = vg.Integrator(DIM * [[0.0, 1.0]])

            integrals = run_vegas(
                batch_f,
                integ,
                NINT_warmup=self.NINT_warmup,
                NEVAL_warmup=self.NEVAL_warmup,
                NINT=self.NINT,
                NEVAL=self.NEVAL,
            )

        #########################
        # Get the int variables and weights
        samples, weights = get_samples(integ, batch_f)

        four_momenta = integrands.get_momenta_from_vegas_samples(
            vsamples=samples, MC_case=self
        )

        # SAVE ALL EVENTS AS A PANDAS DATAFRAME
        df_gen = create_df_from_vegas(four_momenta=four_momenta)

        # add weights to it
        df_gen["w_flux"] = weights["diff_rate"]
        df_gen["w_decay_rate"] = weights["diff_decay_rate"]

        return df_gen


def create_df_from_vegas(four_momenta, sparse=0):
    particles = list(four_momenta.keys())

    if sparse >= 2:  # keep visible, and parent momenta -- Enu to be added later
        particles = [
            x
            for x in particles
            if "target" not in x and "recoils" not in x and "daughter" not in x
        ]
    if sparse == 4:  # keep only visible momenta
        particles = [x for x in particles if "w_decay" not in x]

    columns_index = pd.MultiIndex.from_product([particles, ["0", "1", "2", "3"]])

    df_gen = pd.DataFrame(
        np.hstack([four_momenta[p] for p in particles]), columns=columns_index
    )

    # differential weights
    for column in df_gen:
        if "w_" in str(column):
            df_gen[column, ""] = df_gen[column]
    return df_gen
