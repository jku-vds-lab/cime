import { Box, Typography } from "@material-ui/core";
import { DatasetDrop, LoadingIndicatorDialog, SDFLoader, SDFModifierDialog, useCancellablePromise } from "projection-space-explorer";
import React from "react";
import { UploadedFiles } from "./UploadedFiles";


export function DatasetTabPanel({ onDataSelected }) {
    const [entry, setEntry] = React.useState(null);
    const [openSDFDialog, setOpen] = React.useState(false);
    const [refreshUploadedFiles, setRefreshUploadedFiles] = React.useState(0);

    const { cancellablePromise, cancelPromises } = useCancellablePromise();
    let abort_controller = new AbortController();

    function onModifierDialogClose(modifiers) {
        setOpen(false);
        if (modifiers !== null){
            abort_controller = new AbortController();
            new SDFLoader().resolvePath(entry, onDataSelected, cancellablePromise, modifiers, abort_controller);
        }
    }

    return <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
        <Box paddingLeft={2} paddingTop={2}>
            <Typography variant="subtitle2" gutterBottom>{'Custom Datasets (Drag and Drop)'}</Typography>
        </Box>

        <DatasetDrop onChange={(var1, var2) => {
                onDataSelected(var1, var2);
                setRefreshUploadedFiles(refreshUploadedFiles + 1);
            }} cancellablePromise={cancellablePromise} abort_controller={abort_controller} />



        <Box paddingLeft={2} paddingTop={2}>
            <Typography variant="subtitle2" gutterBottom>{'Predefined Datasets'}</Typography>
        </Box>
        
        <UploadedFiles onChange={(entry)=>{
            setEntry(entry);
            setOpen(true);
        }} refresh={refreshUploadedFiles} />


        <LoadingIndicatorDialog handleClose={() => {cancelPromises();}} area={"global_loading_indicator"}/>
        <SDFModifierDialog openSDFDialog={openSDFDialog} handleClose={onModifierDialogClose}></SDFModifierDialog>
    </div>
}

