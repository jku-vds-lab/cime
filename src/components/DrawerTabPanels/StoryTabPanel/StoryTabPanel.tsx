import { Box, Button, Grid, makeStyles, Typography } from "@material-ui/core";
import {unweighted} from "graphology-shortest-path";
import React = require("react");
import { connect, ConnectedProps } from "react-redux";
import { openStoryEditor } from "../../Ducks/StoryEditorDuck";
import { RootState } from "../../Store/Store";
import { StoryPreview } from "./StoryPreview/StoryPreview";
const Graph = require('graphology');

const mapStateToProps = (state: RootState) => ({
    stories: state.stories,
    storyEditor: state.storyEditor
})

const mapDispatchToProps = dispatch => ({
    openStoryEditor: visible => dispatch(openStoryEditor(visible))
})

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>

type Props = PropsFromRedux & {

}

const useStyles = makeStyles((theme) => ({
    title: {
        margin: theme.spacing(0, 0, 2),
    },
    list: {
        maxHeight: 400,
        overflow: 'auto'
    }
}));


export var StoryTabPanel = connector(({
    stories,
    openStoryEditor,
    storyEditor }: Props) => {
        
    const classes = useStyles()
    

    return <Grid container spacing={0} direction="column" alignItems="stretch">
        <Grid item>
            <Box p={2}>
                <Typography variant="h6" className={classes.title}>Story Books</Typography>
                <StoryPreview></StoryPreview>
            </Box>
        </Grid>
    </Grid>
})