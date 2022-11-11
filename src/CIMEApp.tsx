import * as React from "react";
import {
  PSEContextProvider,
  API,
  Application,
  PluginRegistry,
  createRootReducer,
  BaseConfig,
  FeatureConfig,
  ComponentConfig,
  DEFAULT_UMAP_SETTINGS,
} from "projection-space-explorer";
import { ChemPlugin } from "./Cime/ChemPlugin";
import { DatasetTabPanel } from "./Overrides/DatasetTabPanel";
import { CimeAppBar } from "./Overrides/CimeAppBar";
import { AppState, CIMEReducers } from "./State/Store";
import merge from 'lodash/merge'
import cloneDeep from 'lodash/cloneDeep'
import "./LineUpContext.scss";

export const DEMO = false;

PluginRegistry.getInstance().registerPlugin(new ChemPlugin());


export type CIMEAppProps = {
  config?: BaseConfig;
  features?: FeatureConfig;
  overrideComponents?: ComponentConfig;
  pseRef?: any;
  providePSEContext?: boolean;
}

export const DEFAULT_CIME_APP_CONFIG: CIMEAppProps = {
  config: {
    preselect: {
      initOnMount: false
    }
  },
  features: {
    embeddings: [
      {id:"umap", name:"UMAP", settings: DEFAULT_UMAP_SETTINGS},
    ],
    showSummaryAttributes: false
  },
  overrideComponents: {
    datasetTab: DatasetTabPanel,
    appBar: CimeAppBar,
  }
}

// TODO: @moritz We are currently using the global object directly, ideally we make it passable as prop.
export const CIMEAppContext = new API<AppState>(undefined, createRootReducer(CIMEReducers));

export function CIMEApp({ providePSEContext = true, ...props }: CIMEAppProps) {
  const [merged, setMerged] = React.useState<CIMEAppProps | null>(null);

  React.useEffect(() => {
    setMerged(merge(cloneDeep(DEFAULT_CIME_APP_CONFIG), props))
  }, [props.config, props.features, props.overrideComponents]);

  const component = <div style={{ width: '100%', height: '100%' }}>
    <Application
      ref={props.pseRef}
      config={merged?.config}
      features={merged?.features}
      //@ts-ignore
      overrideComponents={merged?.overrideComponents}
    />
  </div>


  return (providePSEContext ?
    <PSEContextProvider context={CIMEAppContext}>
      {component}
    </PSEContextProvider> : component
  );
}