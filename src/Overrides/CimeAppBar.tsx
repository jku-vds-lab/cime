import * as React from "react";
import { Toolbar } from "@mui/material";
import { AppBar, Typography } from "@mui/material";
import BayerLogo from "../assets/bayer-logo.svg";

export function CimeAppBar() {
  return (
    <AppBar
      variant="outlined"
      elevation={0}
      position="relative"
      color="transparent"
    >
      <Toolbar>
        <a href={"https://www.bayer.com"} target={"_blank"} rel="noreferrer">
          <img style={{ height: 48 }} src={BayerLogo} alt="Powered By Bayer" />
        </a>
        <Typography
          variant="h6"
          style={{ marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" }}
        >
          {"CIME: ChemInformatics Model Explorer"}
        </Typography>
      </Toolbar>
    </AppBar>
  );
}
