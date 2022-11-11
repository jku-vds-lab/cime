/// <reference types="react" />
import { API, BaseConfig, FeatureConfig, ComponentConfig } from "projection-space-explorer";
import { AppState } from "./State/Store";
import "./LineUpContext.scss";
export declare const DEMO = false;
export declare type CIMEAppProps = {
    config?: BaseConfig;
    features?: FeatureConfig;
    overrideComponents?: ComponentConfig;
    pseRef?: any;
    providePSEContext?: boolean;
};
export declare const DEFAULT_CIME_APP_CONFIG: CIMEAppProps;
export declare const CIMEAppContext: API<AppState>;
export declare function CIMEApp({ providePSEContext, ...props }: CIMEAppProps): JSX.Element;
//# sourceMappingURL=CIMEApp.d.ts.map