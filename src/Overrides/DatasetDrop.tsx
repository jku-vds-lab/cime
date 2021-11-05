import { Grid } from "@mui/material";
import { CSVLoader, JSONLoader } from "projection-space-explorer";
import { DragAndDrop } from "projection-space-explorer";
import React from "react";
import { SDFLoader } from "./SDFLoader";
import { SDFModifierDialog } from "./SDFModifierDialog";

export const DatasetDrop = ({
  onChange,
  cancellablePromise,
  abort_controller,
}) => {
  const [entry, setEntry] = React.useState(null);
  const [openSDFDialog, setOpen] = React.useState(false);

  function onModifierDialogClose(modifiers) {
    setOpen(false);
    if (modifiers !== null) {
      abort_controller = new AbortController();
      new SDFLoader().resolveContent(
        entry,
        onChange,
        cancellablePromise,
        modifiers,
        abort_controller
      );
    }
  }

  return (
    <Grid
      container
      item
      alignItems="stretch"
      justifyContent="center"
      direction="column"
      style={{ padding: "16px" }}
    >
      {/**@ts-ignore**/}
      <DragAndDrop
        accept="*"
        handleDrop={(files) => {
          if (files == null || files.length <= 0) {
            return;
          }

          var file = files[0];
          var fileName = file.name as string;

          if (fileName.endsWith("sdf")) {
            setEntry(file);
            setOpen(true);
          } else {
            var reader = new FileReader();
            reader.onload = (event) => {
              var content = event?.target?.result;

              if (fileName.endsWith("json")) {
                new JSONLoader().resolveContent(content, onChange);
              } else {
                new CSVLoader().resolveContent(content, onChange);
              }
            };

            reader.readAsText(file);
          }
        }}
      >
        <div style={{ height: 200 }}></div>
      </DragAndDrop>
      <SDFModifierDialog
        openSDFDialog={openSDFDialog}
        handleClose={onModifierDialogClose}
      ></SDFModifierDialog>
    </Grid>
  );
};
