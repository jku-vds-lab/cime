import { connect, ConnectedProps } from 'react-redux'
import React = require('react')
import { FlexParent } from '../../util/FlexParent'
import { Box, Button } from '@material-ui/core'
import { ProjectionControlCard } from './ProjectionControlCard/ProjectionControlCard'
import { setProjectionOpenAction } from "../../Ducks/ProjectionOpenDuck"
import { setProjectionWorkerAction } from "../../Ducks/ProjectionWorkerDuck"
import { Dataset } from "../../util/Data/Dataset"
import { GenericSettings } from './GenericSettings/GenericSettings'
import { RootState } from '../../Store/Store'
import { setProjectionParamsAction } from '../../Ducks/ProjectionParamsDuck'
import { setProjectionColumns } from '../../Ducks/ProjectionColumnsDuck'
import { TSNEEmbeddingController } from './EmbeddingController/TSNEEmbeddingController'
import { UMAPEmbeddingController } from './EmbeddingController/UMAPEmbeddingController'
import { ClusterTrailSettings } from './ClusterTrailSettings/ClusterTrailSettings'
import { setTrailVisibility } from '../../Ducks/TrailSettingsDuck'
import { ForceAtlas2EmbeddingController } from './EmbeddingController/ForceAtlas2EmbeddingController'

const mapStateToProps = (state: RootState) => ({
    currentAggregation: state.currentAggregation,
    stories: state.stories,
    storyMode: state.storyMode,
    projectionWorker: state.projectionWorker,
    projectionOpen: state.projectionOpen,
    dataset: state.dataset,
    webGLView: state.webGLView,
    projectionParams: state.projectionParams
})

const mapDispatchToProps = dispatch => ({
    setProjectionOpen: value => dispatch(setProjectionOpenAction(value)),
    setProjectionWorker: value => dispatch(setProjectionWorkerAction(value)),
    setProjectionParams: value => dispatch(setProjectionParamsAction(value)),
    setProjectionColumns: value => dispatch(setProjectionColumns(value)),
    setTrailVisibility: visibility => dispatch(setTrailVisibility(visibility))
})

const connector = connect(mapStateToProps, mapDispatchToProps);

type PropsFromRedux = ConnectedProps<typeof connector>

type Props = PropsFromRedux & {
    projectionWorker?: Worker
    projectionOpen?: boolean
    setProjectionOpen?: any
    setProjectionWorker?: any
    dataset?: Dataset
    webGLView?: any
}


export const EmbeddingTabPanel = connector((props: Props) => {
    const [open, setOpen] = React.useState(false)
    const [domainSettings, setDomainSettings] = React.useState('')

    const [controller, setController] = React.useState(null)


    return <FlexParent
        alignItems='stretch'
        flexDirection='column'
        justifyContent=''
    >
        <Box p={1}>
            <ProjectionControlCard
                controller={controller}
                onClose={() => {
                    if (controller) {
                        controller.terminate()
                    }
                    setController(null)
                    props.setTrailVisibility(false)
                }}
                onComputingChanged={(e, newVal) => {
                }} />
        </Box>

        <Box p={1}>
            <Button
                style={{
                    width: '100%'
                }}
                variant="outlined"
                onClick={() => {
                    setDomainSettings('umap')
                    setOpen(true)
                }}>{'UMAP'}</Button>
        </Box>


        <Box p={1}>
            <Button
                style={{
                    width: '100%'
                }}
                variant="outlined"
                onClick={() => {
                    setDomainSettings('tsne')
                    setOpen(true)
                }}>{'t-SNE'}</Button>
        </Box>

        <Box p={1}>
            <Button
                style={{
                    width: '100%'
                }}
                variant="outlined"
                onClick={() => {
                    setDomainSettings('forceatlas2')
                    setOpen(true)
                }}>{'ForceAtlas2'}</Button>
        </Box>


        <GenericSettings
            projectionParams={props.projectionParams}
            domainSettings={domainSettings}
            open={open} onClose={() => setOpen(false)}
            onStart={(params, selection) => {
                setOpen(false)
                props.setProjectionColumns(selection)
                props.setProjectionParams(params)

                switch (domainSettings) {
                    case 'tsne': {
                        let controller = new TSNEEmbeddingController()
                        controller.init(props.dataset, selection, params)
                        controller.stepper = (Y) => {
                            props.dataset.vectors.forEach((vector, i) => {
                                vector.x = Y[i][0]
                                vector.y = Y[i][1]
                            })
                            props.webGLView.current.updateXY()
                            props.webGLView.current.repositionClusters()
                        }

                        setController(controller)
                        break;
                    }

                    case 'umap': {
                        let controller = new UMAPEmbeddingController()
                        let samples = params.useSelection ? props.currentAggregation : props.dataset.vectors

                        controller.init(props.dataset, selection, params, params.useSelection ? samples : undefined)
                        controller.stepper = (Y) => {
                            let source = controller.boundsY(Y)
                            let target = controller.targetBounds



                            samples.forEach((sample, i) => {
                                if (controller.targetBounds) {
                                    sample.x = target.x + ((Y[i][0] - source.x) / source.width) * target.width
                                    sample.y = target.y + ((Y[i][1] - source.y) / source.height) * target.height
                                } else {
                                    sample.x = Y[i][0]
                                    sample.y = Y[i][1]
                                }

                            })


                            props.webGLView.current.updateXY()
                            props.webGLView.current.repositionClusters()
                        }

                        setController(controller)
                        break;
                    }
                    case 'forceatlas2': {
                        let controller = new ForceAtlas2EmbeddingController()
                        controller.init(props.dataset, selection, params)

                        controller.stepper = (Y) => {
                            props.dataset.vectors.forEach((sample, i) => {
                                let idx = controller.nodes[sample.view.duplicateOf].view.meshIndex
                                sample.x = Y[idx].x
                                sample.y = Y[idx].y
                            })
                            props.webGLView.current.updateXY()
                        }

                        setController(controller)
                        break;
                    }

                }
            }}
        ></GenericSettings>


        <ClusterTrailSettings></ClusterTrailSettings>
    </FlexParent>
})