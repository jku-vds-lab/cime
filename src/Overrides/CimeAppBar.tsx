import * as React from "react";
import { Box, IconButton, Menu, MenuItem, Toolbar } from "@mui/material";
import { AppBar, Typography } from "@mui/material";
import BayerLogo from "../assets/bayer-logo.svg";
import { styled, alpha } from '@mui/material/styles';
import MoreIcon from '@mui/icons-material/MoreVert';

export function CimeAppBar() {
  const menuId = 'primary-search-account-menu';
  const [anchorEl, setAnchorEl] = React.useState(null);
  const isMenuOpen = Boolean(anchorEl);

  const handleMenuClose = () => {
    setAnchorEl(null);
  };

  const handleProfileMenuOpen = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const renderMenu = (
    <Menu
      anchorEl={anchorEl}
      anchorOrigin={{
        vertical: 'top',
        horizontal: 'right',
      }}
      id={menuId}
      keepMounted
      transformOrigin={{
        vertical: 'top',
        horizontal: 'right',
      }}
      open={isMenuOpen}
      onClose={handleMenuClose}
    >
      <MenuItem onClick={handleMenuClose}>About</MenuItem>
      <MenuItem onClick={handleMenuClose}>Help</MenuItem>
      <MenuItem onClick={handleMenuClose}>...</MenuItem>
    </Menu>
  );

  return (
    <Box sx={{ flexGrow: 1 }}>
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

          <Box sx={{ flexGrow: 1 }} />

          <IconButton
            size="large"
            edge="end"
            aria-label="account of current user"
            aria-controls={menuId}
            onClick={handleProfileMenuOpen}
            aria-haspopup="true"
            color="inherit"
          >
            <MoreIcon />
          </IconButton>
        </Toolbar>
      </AppBar>
      {renderMenu}
    </Box>
  );
}
