import {
  Button,
  Grid,
  IconButton,
  List,
  ListItem,
  ListItemSecondaryAction,
  ListItemText,
  ListSubheader,
} from "@mui/material";
import DeleteIcon from "@mui/icons-material/Delete";
import { DatasetType, useCancellablePromise } from "projection-space-explorer";
import React from "react";
import RefreshIcon from "@mui/icons-material/Refresh";
import { trackPromise } from "react-promise-tracker";
import { DEMO } from "../CIMEApp";
import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import { LoadingIndicatorView } from "./DatasetTabPanel";

interface IUploadedFile {
  name: string;
  id: number;
}

export const UploadedFiles = ({ onChange, refresh }) => {
  const [files, setFiles] = React.useState<IUploadedFile[]>([]);
  const { cancellablePromise } = useCancellablePromise();

  React.useEffect(() => {
    updateFiles();
  }, [refresh]);

  var handleClick = (entry) => {
    onChange(entry);
  };

  const loading_area = "update_uploaded_files_list";

  function updateFiles() {
    trackPromise(
      cancellablePromise(CIMEBackendFromEnv.getFiles())
        .then((data) => setFiles(data))
        .catch((error) => console.log(error)),
      loading_area
    );
  }

  var handleDelete = (file: IUploadedFile) => {
    cancellablePromise(CIMEBackendFromEnv.deleteFile(file.id))
      .then(() => updateFiles())
      .catch((error) => console.log(error));
  };

  return (
    files && (
      <div>
        <Grid
          item
          style={{ overflowY: "auto", flex: "1 1 auto", maxHeight: "400px" }}
        >
          <List subheader={<li />} style={{ backgroundColor: "white" }}>
            {!DEMO && (
              <ListSubheader>
                Uploaded Files{" "}
                <Button onClick={() => updateFiles()}>
                  <RefreshIcon style={{ fontSize: "1.25rem" }} />
                </Button>
              </ListSubheader>
            )}
            {DEMO && <ListSubheader>Select Dataset</ListSubheader>}
            {files.map((file) => (
              <ListItem
                key={file.id}
                button
                onClick={() => {
                  handleClick({
                    display: file.name,
                    path: file.id,
                    type: DatasetType.Chem,
                    uploaded: true, // indicates that file is already uploaded
                  });
                }}
              >
                <ListItemText primary={file.name}></ListItemText>
                {!DEMO && (
                  <ListItemSecondaryAction
                    onClick={() => {
                      handleDelete(file);
                    }}
                  >
                    <IconButton edge="end" aria-label="delete">
                      <DeleteIcon />
                    </IconButton>
                  </ListItemSecondaryAction>
                )}
              </ListItem>
            ))}
          </List>
          <LoadingIndicatorView area={loading_area} />
        </Grid>
      </div>
    )
  ); //<Divider/>
};
