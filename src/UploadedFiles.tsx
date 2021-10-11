import { Button, Divider, Grid, IconButton, List, ListItem, ListItemSecondaryAction, ListItemText, ListSubheader } from "@material-ui/core"
import DeleteIcon from '@material-ui/icons/Delete';
import { DatasetType, delete_file, get_uploaded_files, LoadingIndicatorView, useCancellablePromise } from "projection-space-explorer";
import React from "react";
import { BiRefresh } from "react-icons/bi";
import { trackPromise } from "react-promise-tracker";
import { DEMO } from "./App";


export const UploadedFiles = ({ onChange, refresh }) => {
    const [files, setFiles] = React.useState<any>(null);
    const { cancellablePromise } = useCancellablePromise();

    React.useEffect(()=>{
        update_files();
    }, [refresh]);

    var handleClick = (entry) => {
        onChange(entry)
    }

    const loading_area = "update_uploaded_files_list";
    
    function update_files(){
        trackPromise(
            cancellablePromise(get_uploaded_files()).then((data: any) => {
                if(data && Object.keys(data).includes("file_list"))
                    setFiles(data["file_list"])
        }).catch(error => console.log(error)), loading_area);
    }
    
    var handleDelete = (file_name) => {
        cancellablePromise(delete_file(file_name)).then((x: any) => {
            if(x && x['deleted'] == "true"){
                const new_files = [...files];
                const index = new_files.indexOf(file_name);
                if (index > -1) {
                    new_files.splice(index, 1);
                }
                setFiles(new_files);
            }
        }).catch(error => console.log(error));
    }
    
    return (files && <div> 
        <Grid item style={{ overflowY: 'auto', flex: '1 1 auto', maxHeight: '400px' }}>

        <List subheader={<li />} style={{ backgroundColor: 'white' }}>
            {!DEMO && <ListSubheader>Uploaded Files <Button onClick={() => update_files()}><BiRefresh/></Button></ListSubheader>}
            {DEMO && <ListSubheader>Select Dataset</ListSubheader>}
            {
            files.map(file_name => (
                <ListItem style={{ maxWidth:'270px', }} key={file_name} button onClick={() => {
                    handleClick({
                        display: file_name,
                        path: file_name,
                        type: DatasetType.Chem,
                        uploaded: true // indicates that file is already uploaded
                    })
                }
                }>
                    <ListItemText primary={file_name}></ListItemText>
                    {!DEMO && <ListItemSecondaryAction onClick={() => {handleDelete(file_name)}}>
                        <IconButton edge="end" aria-label="delete">
                            <DeleteIcon />
                        </IconButton>
                    </ListItemSecondaryAction>}
                </ListItem>
            ))
            }
        </List>
        <LoadingIndicatorView area={loading_area}/>
    </Grid></div>) //<Divider/>
}

