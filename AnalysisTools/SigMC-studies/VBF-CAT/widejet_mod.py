#!/usr/bin/env python3
"""
widejet_modified.py

Modified CMS-style wide-jet algorithm.

This module implements a custom wide-jet reconstruction in which:

  1. Seed 1 is the highest-pT jet in the event.
  2. Seed 2 is the jet closest in phi to the direction opposite seed 1.
  3. Additional jets are clustered into the nearest seed if they lie within
     a configurable wide-jet radius, by default DeltaR < 1.1.
  4. Jets not assigned to either wide jet are retained as "leftover jets".
  5. The leftover jets are used to classify the event as VBF-like if at least
     one leftover-jet pair satisfies |DeltaEta| > 3.

The code assumes that jets are ROOT Lorentz vectors, for example
ROOT.TLorentzVector objects, with methods:

    Pt()
    Eta()
    Phi()
    M()
    Clone()

and support for four-vector addition using += and +.

Author: generated helper implementation
"""

import math


def delta_phi(phi1, phi2):
    """
    Compute the difference in azimuthal angle between two directions.

    Parameters
    ----------
    phi1, phi2 : float
        Azimuthal angles in radians.

    Returns
    -------
    float
        The angular difference phi1 - phi2, folded into the interval
        [-pi, +pi].

    Why this is needed
    ------------------
    Phi is periodic. For example, +pi and -pi correspond to the same
    physical direction. Therefore, a naive subtraction can give a large
    number even when two directions are actually close.
    """

    dphi = phi1 - phi2

    # Fold the angle back into the conventional range.
    while dphi > math.pi:
        dphi -= 2.0 * math.pi

    while dphi <= -math.pi:
        dphi += 2.0 * math.pi

    return dphi


def wrap_phi(phi):
    """
    Fold an azimuthal angle into the interval [-pi, +pi].

    This is useful when we explicitly construct a new direction,
    for example phi + pi, and want to express it in the usual range.
    """

    while phi > math.pi:
        phi -= 2.0 * math.pi

    while phi <= -math.pi:
        phi += 2.0 * math.pi

    return phi


def delta_r(jet1, jet2):
    """
    Compute the standard hadron-collider angular distance DeltaR.

    DeltaR is defined as:

        DeltaR = sqrt( DeltaEta^2 + DeltaPhi^2 )

    Parameters
    ----------
    jet1, jet2 : ROOT Lorentz-vector-like objects
        Objects with Eta() and Phi() methods.

    Returns
    -------
    float
        The DeltaR distance between the two input jets.
    """

    deta = jet1.Eta() - jet2.Eta()
    dphi = delta_phi(jet1.Phi(), jet2.Phi())

    return math.sqrt(deta * deta + dphi * dphi)


def clone_lorentz_vector(vec):
    """
    Clone a ROOT Lorentz vector.

    For ROOT.TLorentzVector, Clone() exists and returns a copy.
    This helper keeps the cloning operation in one place, so that if
    you later use a different Lorentz-vector class, you only need to
    adapt this function.

    Parameters
    ----------
    vec : ROOT Lorentz-vector-like object

    Returns
    -------
    ROOT Lorentz-vector-like object
        A copy of the input four-vector.
    """

    # ROOT.TLorentzVector supports Clone().
    if hasattr(vec, "Clone"):
        return vec.Clone()

    # Fallback for some ROOT.Math Lorentz vectors:
    # their Python bindings often support copy construction via type(vec)(vec).
    try:
        return type(vec)(vec)
    except Exception as exc:
        raise TypeError(
            "Could not clone the input Lorentz vector. "
            "For ROOT.TLorentzVector this should work automatically. "
            "If you are using ROOT.Math LorentzVector, adapt clone_lorentz_vector()."
        ) from exc


def find_seed_indices(jets):
    """
    Find the two seed jets according to the modified prescription.

    Seed 1:
        the highest-pT jet in the event.

    Seed 2:
        the jet closest in phi to the direction opposite seed 1.

    Parameters
    ----------
    jets : list
        List of ROOT Lorentz-vector-like jet objects.

    Returns
    -------
    tuple
        (seed1_index, seed2_index)

    Notes
    -----
    Seed 2 is not necessarily the second-highest-pT jet. It is selected
    geometrically in phi, relative to seed 1.
    """

    if len(jets) < 2:
        raise ValueError("At least two jets are required to define two seeds.")

    # Seed 1 is the highest-pT jet in the event.
    seed1_index = max(range(len(jets)), key=lambda i: jets[i].Pt())
    seed1 = jets[seed1_index]

    # The phi direction opposite to seed 1 is phi(seed1) + pi.
    # We wrap it back into [-pi, +pi].
    opposite_phi = wrap_phi(seed1.Phi() + math.pi)

    # Seed 2 is chosen among all jets except seed 1.
    # We choose the jet with the smallest absolute DeltaPhi from opposite_phi.
    candidate_indices = [i for i in range(len(jets)) if i != seed1_index]

    seed2_index = min(
        candidate_indices,
        key=lambda i: abs(delta_phi(jets[i].Phi(), opposite_phi)),
    )

    return seed1_index, seed2_index


def classify_vbf_from_leftover_jets(leftover_jets, leftover_indices, vbf_delta_eta_min=3.0):
    """
    Classify the event as VBF-like using the leftover jets.

    The event is tagged as VBF-like if there exists at least one pair of
    leftover jets satisfying:

        |eta_i - eta_j| > vbf_delta_eta_min

    Parameters
    ----------
    leftover_jets : list
        Jets not used in the wide-jet reconstruction.

    leftover_indices : list
        Indices of the leftover jets in the original jet collection.

    vbf_delta_eta_min : float
        Minimum pseudorapidity separation required for the VBF tag.

    Returns
    -------
    tuple
        (is_vbf, vbf_pair_indices)

        is_vbf : bool
            True if at least one leftover-jet pair passes the requirement.

        vbf_pair_indices : tuple or None
            Pair of original jet indices corresponding to the first pair found.
            None if no pair satisfies the requirement.
    """

    # We need at least two leftover jets to form a VBF pair.
    if len(leftover_jets) < 2:
        return False, None

    # Try all unique pairs: (0,1), (0,2), ..., (1,2), ...
    for a in range(len(leftover_jets)):
        for b in range(a + 1, len(leftover_jets)):

            deta = abs(leftover_jets[a].Eta() - leftover_jets[b].Eta())

            if deta > vbf_delta_eta_min:
                return True, (leftover_indices[a], leftover_indices[b])

    return False, None


def build_custom_widejets_and_vbf_tag(
    jets,
    widejet_radius=1.1,
    vbf_delta_eta_min=3.0,
):
    """
    Build modified wide jets and classify the event as VBF-like or not.

    Parameters
    ----------
    jets : list
        List of ROOT Lorentz-vector-like reconstructed jets.

        The jets can be in any order. The function internally finds
        the leading-pT jet to define seed 1.

    widejet_radius : float, optional
        Maximum DeltaR from a seed for an additional jet to be absorbed
        into a wide jet. The CMS wide-jet value is usually 1.1.

    vbf_delta_eta_min : float, optional
        Minimum |DeltaEta| required between two leftover jets to tag the
        event as VBF-like. The default is 3.0.

    Returns
    -------
    dict
        Dictionary containing:

        widejet1
            Four-vector of the first wide jet.

        widejet2
            Four-vector of the second wide jet.

        seed1_index
            Index of seed 1 in the original jet list.

        seed2_index
            Index of seed 2 in the original jet list.

        assigned_to_widejet1_indices
            Original indices of jets included in widejet1.

        assigned_to_widejet2_indices
            Original indices of jets included in widejet2.

        leftover_jets
            List of jets not used in either wide jet.

        leftover_indices
            Original indices of the leftover jets.

        is_vbf
            Boolean VBF classification.

        vbf_pair_indices
            Original jet indices of the first leftover pair satisfying
            |DeltaEta| > vbf_delta_eta_min, or None.
    """

    # If fewer than two jets are present, the wide-jet system cannot be built.
    # We return a complete dictionary anyway, so downstream code can handle
    # this event without crashing.
    if len(jets) < 2:
        return {
            "widejet1": None,
            "widejet2": None,
            "seed1_index": None,
            "seed2_index": None,
            "assigned_to_widejet1_indices": [],
            "assigned_to_widejet2_indices": [],
            "leftover_jets": list(jets),
            "leftover_indices": list(range(len(jets))),
            "is_vbf": False,
            "vbf_pair_indices": None,
        }

    # Find the two seed jets according to the modified seed prescription.
    seed1_index, seed2_index = find_seed_indices(jets)

    seed1 = jets[seed1_index]
    seed2 = jets[seed2_index]

    # Initialize the wide jets as copies of the two seed four-vectors.
    # We clone them to avoid modifying the original jet objects.
    widejet1 = clone_lorentz_vector(seed1)
    widejet2 = clone_lorentz_vector(seed2)

    # Keep explicit bookkeeping of which original jets enter each wide jet.
    # The seed jets are part of the corresponding wide jets by definition.
    assigned_to_widejet1_indices = [seed1_index]
    assigned_to_widejet2_indices = [seed2_index]

    # Jets not assigned to either wide jet are retained here.
    leftover_jets = []
    leftover_indices = []

    # Loop over all reconstructed jets in the event.
    for i, jet in enumerate(jets):

        # Skip the seeds because they are already included.
        if i == seed1_index or i == seed2_index:
            continue

        # Compute the distance of this jet from each original seed axis.
        #
        # Important: we compare to the original seed directions, not to the
        # progressively updated wide-jet four-vectors. This keeps the axes
        # stable and close to the CMS wide-jet construction.
        dR1 = delta_r(jet, seed1)
        dR2 = delta_r(jet, seed2)

        # Case 1:
        # The jet is inside the wide-jet radius around seed 1 and is closer
        # to seed 1 than to seed 2. Add it to widejet1.
        if dR1 < widejet_radius and dR1 < dR2:
            widejet1 += jet
            assigned_to_widejet1_indices.append(i)

        # Case 2:
        # The jet is inside the wide-jet radius around seed 2 and is closer
        # to seed 2 than to seed 1. Add it to widejet2.
        elif dR2 < widejet_radius and dR2 < dR1:
            widejet2 += jet
            assigned_to_widejet2_indices.append(i)

        # Case 3:
        # The jet is outside both wide-jet cones, or exactly tied in distance.
        # We do not use it in the wide-jet system and retain it as a leftover.
        else:
            leftover_jets.append(jet)
            leftover_indices.append(i)

    # Use the leftover jets to classify the event as VBF-like or not.
    is_vbf, vbf_pair_indices = classify_vbf_from_leftover_jets(
        leftover_jets=leftover_jets,
        leftover_indices=leftover_indices,
        vbf_delta_eta_min=vbf_delta_eta_min,
    )

    return {
        "widejet1": widejet1,
        "widejet2": widejet2,
        "seed1_index": seed1_index,
        "seed2_index": seed2_index,
        "assigned_to_widejet1_indices": assigned_to_widejet1_indices,
        "assigned_to_widejet2_indices": assigned_to_widejet2_indices,
        "leftover_jets": leftover_jets,
        "leftover_indices": leftover_indices,
        "is_vbf": is_vbf,
        "vbf_pair_indices": vbf_pair_indices,
    }


def print_result_summary(result):
    """
    Print a compact summary of the wide-jet reconstruction result.
    """

    print("=== Modified wide-jet result ===")

    print(f"Seed 1 index: {result['seed1_index']}")
    print(f"Seed 2 index: {result['seed2_index']}")

    print(f"Jets in widejet1: {result['assigned_to_widejet1_indices']}")
    print(f"Jets in widejet2: {result['assigned_to_widejet2_indices']}")

    print(f"Leftover jet indices: {result['leftover_indices']}")

    print(f"VBF-like event: {result['is_vbf']}")
    print(f"VBF pair indices: {result['vbf_pair_indices']}")

    if result["widejet1"] is not None and result["widejet2"] is not None:
        mjj = (result["widejet1"] + result["widejet2"]).M()
        print(f"Wide-jet dijet mass: {mjj:.3f}")


if __name__ == "__main__":
    """
    Example usage.

    This block is executed only if you run:

        python widejet_modified.py

    In a real analysis, you would usually import the function from this file:

        from widejet_modified import build_custom_widejets_and_vbf_tag

    and then call it inside your event loop.
    """

    try:
        import ROOT
    except ImportError:
        raise ImportError(
            "This example requires PyROOT. "
            "Run it in a ROOT-enabled environment, or import the functions "
            "from this file inside your existing CMSSW/ROOT analysis."
        )

    # ------------------------------------------------------------------
    # Create a toy list of ROOT.TLorentzVector jets.
    #
    # In your real analysis, replace this block with the jets read from
    # your event tree.
    #
    # For example, your real code may look like:
    #
    #     jets = []
    #     for i in range(nJet):
    #         v = ROOT.TLorentzVector()
    #         v.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i])
    #         jets.append(v)
    #
    # ------------------------------------------------------------------

    jets = []

    # Leading jet: this will become seed 1.
    j0 = ROOT.TLorentzVector()
    j0.SetPtEtaPhiM(500.0, 0.3, 0.2, 50.0)
    jets.append(j0)

    # Jet approximately opposite in phi to seed 1: this will become seed 2.
    j1 = ROOT.TLorentzVector()
    j1.SetPtEtaPhiM(420.0, -0.2, wrap_phi(0.2 + math.pi + 0.1), 45.0)
    jets.append(j1)

    # Nearby radiation around seed 1: should be added to widejet1.
    j2 = ROOT.TLorentzVector()
    j2.SetPtEtaPhiM(80.0, 0.5, 0.45, 10.0)
    jets.append(j2)

    # Nearby radiation around seed 2: should be added to widejet2.
    j3 = ROOT.TLorentzVector()
    j3.SetPtEtaPhiM(60.0, -0.4, wrap_phi(0.2 + math.pi - 0.2), 8.0)
    jets.append(j3)

    # Leftover forward jet.
    j4 = ROOT.TLorentzVector()
    j4.SetPtEtaPhiM(70.0, 3.2, 1.5, 12.0)
    jets.append(j4)

    # Another leftover forward/backward jet.
    # Together with j4, this gives |DeltaEta| > 3 and should tag VBF.
    j5 = ROOT.TLorentzVector()
    j5.SetPtEtaPhiM(65.0, -1.0, -1.0, 10.0)
    jets.append(j5)

    # ------------------------------------------------------------------
    # Run the modified wide-jet reconstruction.
    # ------------------------------------------------------------------

    result = build_custom_widejets_and_vbf_tag(
        jets,
        widejet_radius=1.1,
        vbf_delta_eta_min=3.0,
    )
    
    print(result)

    # ------------------------------------------------------------------
    # Print the result.
    # ------------------------------------------------------------------

    print_result_summary(result)
