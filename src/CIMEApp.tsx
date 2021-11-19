import * as React from "react";
import { useState } from "react";
import {
  PSEContextProvider,
  API,
  Application,
  PluginRegistry,
  createRootReducer,
  PSEIcons,
  BaseConfig,
  FeatureConfig,
  ComponentConfig,
  DatasetType,
  setDatasetEntriesAction,
} from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { LineUpContext } from "./LineUpContext";
import { LineUpTabPanel } from "./Overrides/LineUpTabPanel";
import { AppState, CIMEReducers } from "./State/Store";
import merge from 'lodash/merge'
import clone from "lodash/clone";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ChemPlugin());


export type CIMEAppProps = {
  config?: BaseConfig
  features?: FeatureConfig
  overrideComponents?: ComponentConfig
}

export const DEFAULT_CIME_APP_CONFIG = {
  config: {
    preselect: {
      initOnMount: false
    }
  },
  features: {
    disableEmbeddings: {
      tsne: true,
      forceatlas: true,
    },
  },
  overrideComponents: {
    datasetTab: DatasetTabPanel,
    appBar: CimeAppBar,
    tabs: [
      {
        name: "lineup",
        tab: LineUpTabPanel,
        title: "LineUp Integration",
        description: "Settings for LineUp Integration",
        icon: PSEIcons.PseLineup,
      },
    ],
    detailViews: [
      {
        name: "lineup",
        view: LineUpContext,
      },
    ],
  }
}


const context = new API<AppState>(undefined, createRootReducer(CIMEReducers))


export function CIMEApp(props: CIMEAppProps) {
  const [merged, setMerged] = React.useState(merge(clone(DEFAULT_CIME_APP_CONFIG), props))

  React.useEffect(() => {
    const cimeAppConfig = clone(DEFAULT_CIME_APP_CONFIG)

    setMerged(merge(cimeAppConfig, props))
  }, [props.config, props.features, props.overrideComponents])

  return (
    <PSEContextProvider context={context}>
      <div style={{ width: '100%', height: '100%' }}>
        <Application
          config={merged.config}
          features={merged.features}
          //@ts-ignore
          overrideComponents={merged.overrideComponents}
        />
      </div>
    </PSEContextProvider>
  );
}