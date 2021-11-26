import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  Typography,
} from "@mui/material";
import { useCancellablePromise, RootActions } from "projection-space-explorer";
import React from "react";
import { usePromiseTracker } from "react-promise-tracker";
import { DatasetDrop } from "./DatasetDrop";
import { SDFLoader } from "./SDFLoader";
import { SDFModifierDialog } from "./SDFModifierDialog";
import { UploadedFiles } from "./UploadedFiles";
import Loader from "react-loader-spinner";
import { connect, ConnectedProps } from "react-redux";

export const LoadingIndicatorView = (props) => {
  const { promiseInProgress } = usePromiseTracker({ area: props.area });
  return (
    promiseInProgress && (
      <div
        style={{
          width: "100%",
          height: "100",
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
        }}
      >
        <Loader type="ThreeDots" color="#2BAD60" height="100" width="100" />
      </div>
    )
  );
};

export const LoadingIndicatorDialog = (props) => {
  const { promiseInProgress } = usePromiseTracker({ area: props.area });

  return (
    <Dialog maxWidth="lg" open={promiseInProgress}>
      {" "}
      {/*onClose={props.handleClose} */}
      <DialogContent>
        <LoadingIndicatorView area={props.area} />
      </DialogContent>
      <DialogActions>
        <Button onClick={props.handleClose}>Cancel</Button>
      </DialogActions>
    </Dialog>
  );
};


const mapStateToProps = (state: any) => ({
})

const mapDispatchToProps = dispatch => ({
  setDataset: value => dispatch(RootActions.loadDataset(value))
})

const connector = connect(mapStateToProps, mapDispatchToProps);

type Props = ConnectedProps<typeof connector>

export const DatasetTabPanel = connector(({ setDataset }: Props) => {
  const [entry, setEntry] = React.useState(null);
  const [openSDFDialog, setOpen] = React.useState(false);
  const [refreshUploadedFiles, setRefreshUploadedFiles] = React.useState(0);

  const { cancellablePromise, cancelPromises } = useCancellablePromise();
  let abort_controller = new AbortController();

  function onModifierDialogClose(modifiers) {
    //console.log(entry);
    setOpen(false);
    if (modifiers !== null) {
      abort_controller = new AbortController();
      new SDFLoader().resolvePath(
        entry,
        (dataset) => {
          //console.log("Dataset!!!", dataset);

          setDataset(dataset);
        },
        cancellablePromise,
        modifiers,
        abort_controller
      );
    }
  }

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <Box paddingLeft={2} paddingTop={2}>
        <Typography variant="subtitle2" gutterBottom>
          {"Custom Datasets (Drag and Drop)"}
        </Typography>
      </Box>

      <DatasetDrop
        onChange={(dataset) => {
          //console.log("Dataaset", dataset);
          setDataset(dataset);
          setRefreshUploadedFiles(refreshUploadedFiles + 1);
        }}
        cancellablePromise={cancellablePromise}
        abort_controller={abort_controller}
      />

      <Box paddingLeft={2} paddingTop={2}>
        <Typography variant="subtitle2" gutterBottom>
          {"Predefined Datasets"}
        </Typography>
      </Box>

      <UploadedFiles
        onChange={(entry) => {
          setEntry(entry);
          setOpen(true);
        }}
        refresh={refreshUploadedFiles}
      />

      <LoadingIndicatorDialog
        handleClose={() => {
          cancelPromises();
        }}
        area={"global_loading_indicator"}
      />

      <SDFModifierDialog
        openSDFDialog={openSDFDialog}
        handleClose={onModifierDialogClose}
      ></SDFModifierDialog>
    </div>
  );
})