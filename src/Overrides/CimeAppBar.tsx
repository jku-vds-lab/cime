import { Toolbar } from "@mui/material";
import { AppBar, Typography } from "@mui/material";
import { PseAppBar } from "projection-space-explorer";

function JJ({ children }) {
    return <AppBar variant="outlined" position="relative" color="transparent">
        <Toolbar>
            {children}
        </Toolbar>
    </AppBar>
}

export function CimeAppBar() {
    return <JJ>
        <a href={"https://www.bayer.com"} target={"_blank"}><img style={{ height: 48 }} src={"textures/bayer-logo.svg"} alt="Powered By Bayer" /></a>
        <Typography variant="h6" style={{ marginLeft: 48, color: "rgba(0, 0, 0, 0.54)" }}>
            {"CIME: Chem-Informatics Model Explorer"}
        </Typography>
    </JJ>
}