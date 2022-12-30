def find_highres_sphere(self) -> Tuple[np.ndarray, float]:
    """
    Determine the centre and radius of high-res sphere from Velociraptor.

    The selection is made based on the location of the halo in the
    catalogue, optionally after sorting them by M200c or M500c. This
    is determined by the value of `self.params['sort_rule']`.

    This function is only executed by the root rank.

    Returns
    -------
    centre : ndarray(float)
        3-element array holding the centre of the high-res region.
    radius : float
        The target radius of the high-res region, including any requested
        padding.

    """
    # Make sure that we are on the root rank if over MPI
    if comm_rank != 0:
        raise ValueError(f"find_highres_sphere() called on MPI rank {comm_rank}!")

    # Look up the target halo index in Velociraptor catalogue
    vr_index = self.find_halo_index()

    with h5py.File(self.params["vr_file"], "r") as vr_file:

        # First, determine the radius of the high-res region
        r_r200 = 0
        r_r500 = 0
        r200 = vr_file["R_200crit"][vr_index]
        if self.params["highres_radius_r200"] > 0:
            r_r200 = r200 * self.params["highres_radius_r200"]
        try:
            r500 = vr_file["SO_R_500_rhocrit"][vr_index]
            r_r500 = r500 * self.params["highres_radius_r500"]
        except KeyError:
            r500 = None
            if self.params["highres_radius_r500"] > 0:
                warn(
                    "Could not load R500c, ignoring request for "
                    f"minimum high-res radius of "
                    f"{self.params['highres_radius_r500']} r_500.",
                    RuntimeWarning,
                )

        r_highres = max(r_r200, r_r500)
        if r_highres <= 0:
            raise ValueError(f"Invalid radius of high-res region ({r_highres})")

        # If enabled, add a fixed "padding" radius to the high-res sphere
        if self.params["highres_radius_padding"] > 0:
            r_highres += self.params["highres_radius_padding"]

        # If enabled, expand radius to requested minimum.
        if self.params["highres_radius_min"] > 0:
            r_highres = max(r_highres, self.params["highres_radius_min"])

        # Load halo centre
        names = ["X", "Y", "Z"]
        centre = np.zeros(3)
        for icoord, prefix in enumerate(names):
            centre[icoord] = vr_file[f"{prefix}cminpot"][vr_index]

    r500_str = "" if r500 is None else f"{r500:.4f}"
    m200_str = "" if getattr(self, "m200crit", None) is None else f"{self.m200crit:.4f}"
    m500_str = "" if getattr(self, "m500crit", None) is None else f"{self.m500crit:.4f}"
    print(
        "Velociraptor search results:\n"
        f"- Run name: {self.params['fname']}\t"
        f"GroupNumber: {self.params['group_number']}\t"
        f"VR index: {vr_index}\n"
        f"- Centre: {centre[0]:.3f} / {centre[1]:.3f} / {centre[2]:.3f} "
        f"- High-res radius: {r_highres:.4f}\n"
        f"- R_200crit: {r200:.4f}\n"
        f"- R_500crit: {r500_str}\n"
        f"- M_200crit: {m200_str}\n"
        f"- M_500crit: {m500_str}\n"
    )

    return centre, r_highres


def find_halo_index(self) -> int:
    """
    Find the index of the desired target halo.

    This function looks up the desired (field) halo if the selection
    is specified in terms of a position in the mass-ordered list.
    It should only ever be run on the root node, an error is raised if
    this is not the case.

    If the parameter file instructs to sort by M500c, but this is not
    recorded in the Velociraptor catalogue, an error is raised.

    Parameters
    ----------
    None

    Returns
    -------
    halo_index : int
        The catalogue index of the target halo.

    """
    if comm_rank != 0:
        raise ValueError("find_halo_index() called on rank {comm_rank}!")

    # If the parameter file already specified the VR index, we are done
    if self.params["sort_rule"].lower() == "none":
        return self.params["group_number"]

    # ... otherwise, need to load the desired mass type of all (central)
    # VR haloes, sort them, and find the entry we want
    with h5py.File(self.params["vr_file"], "r") as vr_file:
        structType = vr_file["/Structuretype"][:]
        field_halos = np.where(structType == 10)[0]

        sort_rule = self.params["sort_rule"].lower()
        if sort_rule == "m200crit":
            m_halo = vr_file["/Mass_200crit"][field_halos]
        elif sort_rule == "m500crit":
            # If M500 cannot be loaded, an error will be raised
            m_halo = vr_file["/SO_Mass_500_rhocrit"][field_halos]
        else:
            raise ValueError("Unknown sorting rule '{sort_rule}'!")

    # Sort groups by specified mass, in descending order
    sort_key = np.argsort(-m_halo)
    halo_index = sort_key[self.params["group_number"]]

    # Store mass of target halo used for sorting, for later use
    setattr(self, sort_rule, m_halo[halo_index])
    return halo_index
